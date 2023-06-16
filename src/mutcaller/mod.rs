/**
this module handles the UNALIGNED functions in mutcaller
**/

extern crate simplelog;
extern crate csv;
extern crate clap;
extern crate bam;
extern crate serde;
extern crate fastq;
extern crate itertools;


use clap::{App, load_yaml};
use std::str;
use std::error::Error;
use serde::Deserialize;
use std::fmt; 
use csv::ReaderBuilder;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use itertools::Itertools;
use flate2::GzBuilder;
use flate2::Compression;
use simple_log::info;
use std::path::Path;
use fastq::parse_path;
use fastq::each_zipped;
use fastq::RefRecord;
use crate::mutcaller::fastq::Record;
use flate2::{read};
use std::process::{Command, Stdio};
use std::ffi::OsStr;
use std::fs;
use flate2::read::MultiGzDecoder;
use std::time::{Instant};
#[cfg(not(feature = "paris"))]
use log::*;
use simplelog::{Config, WriteLogger, CombinedLogger, LevelFilter};
use crate::countbam::get_current_working_dir;
use crate::countbam::cleanup;
use crate::vcf::{guess_vcf, guess_compression, read_vcf_compressed, read_vcf_uncompressed};

#[derive(Deserialize)]
#[derive(Debug)]
#[derive(Clone)]
pub struct Variant {
    pub seq: String,
    pub start: String,
    pub ref_nt: String,
    pub query_nt: String,
    pub name: String,
    pub class: Option<VariantClass>,
}

#[derive(Deserialize)]
#[derive(Debug)]
#[derive(Clone)]
#[derive(PartialEq)]
pub enum VariantClass {
    SNV,
    MNV,
    Insertion,
    Deletion
}

// Implement `Display` for `Variant`.
impl fmt::Display for Variant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use `self.number` to refer to each positional data point.
        write!(f, "seq: {} start: {} ref_nt: {} query_nt: {} name: {} class: {:?}", self.seq, self.start, self.ref_nt, self.query_nt, self.name, self.class.as_ref().unwrap())
    }
}


#[derive(Debug)]
pub struct Paramsm <'a> {
    pub fastq1: String,
    pub fastq2: String,
    pub genome: String,
    pub bcs: String,
    pub umi_len: usize,
    pub cb_len: usize,
    pub threads: usize,
    pub aligner: String,
    pub aligner_loc: String,
    pub aligner_args: Box<[&'a str]>,
    pub variants: String,
    pub read_len: usize,
    pub output_path: Box<Path>,
    pub keep: bool,
    pub verbose: bool,
    pub vcf: bool,
    pub qual: f64,
}



fn load_params() -> Paramsm <'static>{
    let yaml = load_yaml!("../cli.yml");
    let matches = App::from_yaml(yaml).get_matches();
    let params = matches.subcommand_matches("UNALIGNED").unwrap();
    let fastq1 = params.value_of("fastq1").unwrap();
    let fastq2 = params.value_of("fastq2").unwrap();
    let output = params.value_of("output").unwrap_or("out");
    let genome = params.value_of("genome").unwrap();
    let bcs = params.value_of("barcodes_file").unwrap_or("/Users/sfurlan/develop/mutCaller/data/737K-august-2016.txt.gz");
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap();
    let umi_len = params.value_of("umi_len").unwrap_or("10");
    let umi_len = umi_len.to_string().parse::<u8>().unwrap();
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<u8>().unwrap();
    let read_len = params.value_of("read_len").unwrap_or("90");
    let read_len = read_len.to_string().parse::<usize>().unwrap();
    let aligner = params.value_of("aligner").unwrap_or("minimap2");
    let mut aligner_loc = aligner;
    let aligner_args = ["--sr", "--splice"];
    let qual = params.value_of("qual").unwrap_or("95.0");
    let qual = qual.to_string().parse::<f64>().unwrap();
    if params.is_present("aligner_loc") {
        aligner_loc = params.value_of("aligner_loc").unwrap();
    }
    // if params.is_present("aligner_args") {
    //     aligner_args = params.value_of("aligner_args").unwrap();
    // }
    let variantstring = params.value_of("variants").unwrap();
    let mut _verbose = true;
    if params.is_present("quiet") {
            _verbose = false
    };
    let mut keep = false;
    if params.is_present("keep_files") {
            keep = true
    };
    let outpath = Path::new(&output);
    let params = Paramsm{
        fastq1: fastq1.to_string(),
        fastq2: fastq2.to_string(),
        genome: genome.to_string(),
        output_path: outpath.into(), 
        bcs: bcs.to_string(),
        threads: threads as usize,
        umi_len: umi_len as usize,
        cb_len: cb_len as usize,
        aligner: aligner.to_string(),
        aligner_loc: aligner_loc.to_string(),
        aligner_args:  Box::new(aligner_args),
        variants: variantstring.to_string(),
        read_len: read_len as usize,
        keep: keep,
        verbose: _verbose,
        vcf: false,
        qual: qual
    };
    return check_params(params).unwrap()
}



pub fn check_params(params: Paramsm) -> Result<Paramsm, Box<dyn Error>>{
    let _cu = cleanup(&params.output_path.join("mutcaller.log"), false);
    let log_file_path = &params.output_path.join("mutcaller.log");
    let log_file = log_file_path.to_str().unwrap();
    info!("\n\n\tStarting!\n");
    let wdpb= get_current_working_dir().unwrap();
    let wdir = wdpb.to_str().unwrap();
    info!("\n\n\tCurrent working directory: '{}'\n", wdir);
    if params.verbose {
        eprintln!("\n\nCurrent working directory: '{}'", wdir);
    }
    if params.output_path.is_relative(){
        let a1 = Path::new(wdir).join(&params.output_path);
        let abs_outpath = a1.to_str().unwrap();
        if params.output_path.exists() {
                info!("\n\n\tFound existing output directory: '{}'\n", &abs_outpath);
                warn!("\n\n\t{}", "Existing data in this folder could be lost!!!\n");
               if params.verbose {
                    eprintln!("Found existing output directory: '{}'", &abs_outpath);
                    eprintln!("\t{}", "Existing data in this folder could be lost!!!");
                }
        } else {
            info!("\n\n\tCreating output directory: '{}'\n", &abs_outpath);
            if params.verbose {
                eprintln!("Creating output directory: '{}'", &abs_outpath);
            }
            fs::create_dir(&params.output_path)?;
        }
    } else {
        let abs_outpath = &params.output_path.to_str().unwrap();
        if params.output_path.exists() {
            info!("\n\n\tFound existing output directory: '{}'\n", &abs_outpath);
            warn!("\n\n\t{}", "Existing data in this folder could be lost!!!\n");
            if params.verbose {
                eprintln!("Found existing output directory: '{}'", &abs_outpath);
                eprintln!("\t{}", "Existing data in this folder could be lost!!!");
            }
        } else {
            info!("\n\n\tCreating output directory: '{}'\n", &abs_outpath);
            if params.verbose {
                eprintln!("Creating output directory: '{}'", &abs_outpath);
            }
            fs::create_dir(&params.output_path)?;
        }
    }
    CombinedLogger::init(vec![
        #[cfg(not(feature = "termcolor"))]
        WriteLogger::new(
            LevelFilter::Info,
            Config::default(),
            File::create(log_file).unwrap(),
        ),
    ])
    .unwrap();

    // check if variants file is vcf
    let is_vcf = guess_vcf(&params.variants);
    Ok(Paramsm{
            fastq1: params.fastq1,
            fastq2: params.fastq2,
            genome: params.genome,
            output_path: params.output_path, 
            bcs: params.bcs,
            threads: params.threads,
            umi_len: params.umi_len,
            cb_len: params.cb_len,
            aligner: params.aligner,
            aligner_loc: params.aligner_loc,
            aligner_args: params.aligner_args,
            variants: params.variants,
            read_len: params.read_len,
            keep: params.keep,
            verbose: params.verbose,
            vcf: is_vcf.unwrap(),
            qual: params.qual
    })
}


pub fn read_csv(params: Option<&Paramsm>, file_i: Option<String>, verbose_i: Option<bool>) -> Result<Vec<Variant>, Box<dyn Error>> {
    // Build the CSV reader and iterate over each record.
//     let data = "\
// seq\tstart\tref_nt\tquery_nt\tname
// chr12\t112450407\tA\tG\tPTPN11_227A>G
// chr2\t208248389\tG\tA\tIDH1_132G>A
// chr17\t7674220\tC\tT\tTP53_248C>T";
    let file = {
        if params.is_some(){
            &params.unwrap().variants
        } else {
            file_i.as_ref().unwrap()
        }
    };

    let verbose = {
        if params.is_some(){
            &params.unwrap().verbose
        } else {
            verbose_i.as_ref().unwrap()
        }
    };

    #[derive(Deserialize)]
    struct PreVariant {
        seq: String,
        start: String,
        ref_nt: String,
        query_nt: String,
        name: String
    }
    let path = Path::new(&file);
    let compression = {
        if path.extension().is_some() && path.extension().unwrap() == "gz" {
                true
            }else{
                false
        }
    };
    if compression {
        if *verbose {
            eprintln!("Opening variants file: {}\n", file.to_string());
        }
        info!("Opening variants file: {}\n", file.to_string());
        // let file = File::open(file.to_string()).unwrap();
        let reader = BufReader::new(MultiGzDecoder::new(File::open(file.to_string()).unwrap()));
        let mut rdr = ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_reader(reader);
        let mut csvdata = Vec::new();
        for result in rdr.deserialize() {
            let prevariant: PreVariant = result?;
            csvdata.push(Variant{
                seq: prevariant.seq,
                start: prevariant.start,
                ref_nt: prevariant.ref_nt,
                query_nt: prevariant.query_nt,
                name: prevariant.name,
                class: None,
            });
        }
        Ok(csvdata)
    } else {
        if *verbose {
            eprintln!("Opening variants file: {}\n", file.to_string());
        }
        info!("Opening variants file: {}\n", file.to_string());
        let file = File::open(file.to_string()).unwrap();
        let reader = BufReader::new(file);
        let mut rdr = ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_reader(reader);
        let mut csvdata = Vec::new();
        for result in rdr.deserialize() {
            let prevariant: PreVariant = result?;
            csvdata.push(Variant{
                seq: prevariant.seq,
                start: prevariant.start,
                ref_nt: prevariant.ref_nt,
                query_nt: prevariant.query_nt,
                name: prevariant.name,
                class: None,
            });
        }
        Ok(csvdata)
    }

}




pub fn mutcaller_run() {
    let start = Instant::now();

    let params = load_params();
    info!("\n\n\tParsing Parameters!\n");
    if params.verbose {
        eprintln!("\n\nParsing Parameters!\n");
    }
    info!("\n\n\tChecking programs and parsing variants!\n");
    if params.verbose {
        eprintln!("\n\nChecking programs and parsing variants!\n");
    }
    let _prog_test_res = test_progs(&params);
    let csvdata = {
        if params.vcf {
            let is_compressed = guess_compression(&params.variants);
            if is_compressed.unwrap() {
                read_vcf_compressed(&params.variants, &params.qual, &params.verbose)
            } else {
                read_vcf_uncompressed(&params.variants, &params.qual, &params.verbose)
            }
        } else {
            Ok(read_csv(Some(&params), None, None).unwrap())
        }
    };
    info!("\n\n\tRunning with {} thread(s)!\n", &params.threads);
    if params.verbose {
        eprintln!("\n\nRunning with {} thread(s)!\n", &params.threads);
        // eprintln!("Paramsm: {:?} ", &params);
    }
    let _fqr = fastq(&params);
    info!("done!");
    let _ar = align(&params);
    let mut count_vec = Vec::new();
    for variant in csvdata.unwrap() {
        let classified_variant = classify_variant(&variant);
        if params.verbose {
            eprintln!("\nCorrectly parsed variant: {}", classified_variant.as_ref().unwrap());
        }
        info!("\n\n\tCorrectly parsed variant: {}\n", classified_variant.as_ref().unwrap());
        let data = count_variants_helper(&params, classified_variant.unwrap());
        if data.is_some() {
            count_vec.push(data.unwrap());
        }
        // for svariant in variant{
        //     count_vec.push(count_variants(&params, svariant.clone()));
        // }
    }
    let _none = writer_fn(count_vec, &params);
    let duration = start.elapsed();
    info!("\n\n\tDone!!\n");
    info!("\n\n\tTime elapsed is: {:?}\n", duration);
    if params.verbose{
        eprintln!("\n\nDone!!");
        eprintln!("\n\nTime elapsed is: {:?}", duration);
    }
}

fn classify_variant (variant: &Variant) -> Result<Variant, Box<dyn Error>> {
    let mut classified_variant = variant.clone();
    if classified_variant.ref_nt.len() == 1 && classified_variant.query_nt.len() == 1{
        classified_variant.class = Some(VariantClass::SNV);
        return Ok(classified_variant)
    }
    if classified_variant.ref_nt.len() > 1 && classified_variant.query_nt.len() > 1{
        classified_variant.class = Some(VariantClass::MNV);
        return Ok(classified_variant)
    }
    if classified_variant.ref_nt.len() == 1 && classified_variant.query_nt.len() >= 1{
        classified_variant.class = Some(VariantClass::Insertion);
        return Ok(classified_variant)
    }
    classified_variant.class = Some(VariantClass::Deletion);
    return Ok(classified_variant)
    // return Err(e)
    // return Ok(classified_variant)
}

// minimap2 --MD -a $fa -t 8 mutcaller_R1.fq.gz -o Aligned.mm2.sam
// samtools sort -@ 8 -o Aligned.mm2.sorted.sam Aligned.mm2.sam
// samtools view -b -@ 8 -o Aligned.mm2.sorted.sam Aligned.mm2.sorted.bam
// samtools index -@ 8 Aligned.mm2.sorted.bam

fn test_progs (params: &Paramsm) -> Result<(), Box<dyn Error>>{
    let _output = Command::new(&params.aligner_loc)
                    .arg("-h")
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                     .output()
                     .expect("\n\n*******Failed to execute aligner*******\n\n");
    let _output = Command::new("samtools")
                    .arg("-h")
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                     .output()
                     .expect("\n\n*******Failed to execute samtools*******\n\n");
    // eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    Ok(())
}


fn align (params: &Paramsm)-> Result<(), Box<dyn Error>> {
    let align_sam = &params.output_path.join("Aligned.sam").clone().to_owned();
    let align_sorted_sam = &params.output_path.join("Aligned.sorted.sam").clone().to_owned();
    let align_sorted_bam = &params.output_path.join("Aligned.sortedByCoord.out.bam").clone().to_owned();
    let fastq = &params.output_path.join("mutcaller_R1.fq.gz").clone().to_owned();
    let outfolder = &params.output_path.join("").clone().to_owned();
    if params.aligner == "minimap2" {
        if params.verbose {
        eprintln!("{}", "Aligning reads using minimap2");
        }
        info!("{}", "Aligning reads using minimap2");
        let output = Command::new(&params.aligner_loc)
                        .args(["--MD", "-Y"])
                        .args(&*params.aligner_args)
                        .arg("-a")
                        .arg(params.genome.to_string())
                        .arg("-t")
                        .arg(params.threads.to_string())
                        .arg(fastq.to_str().unwrap())
                        .arg("-o")
                        .arg(align_sam.to_str().unwrap())
                        .stderr(Stdio::piped())
                        .stdout(Stdio::piped())
                        .output()
                        .expect("\n\n*******Failed to execute minimap2*******\n\n");

        if params.verbose {
            eprintln!("{}", String::from_utf8_lossy(&output.stderr));
            eprintln!("{}", "Minimap2 complete; Running samtools sort");
        }
        info!("{}", String::from_utf8_lossy(&output.stderr));
        info!("{}", "Minimap2 complete; Running samtools sort");
        let output = Command::new("samtools")
                    .arg("sort")
                    .arg("-@")
                    .arg(params.threads.to_string())
                    .arg("-o")
                    .arg(align_sorted_sam.to_str().unwrap())
                    .arg(align_sam.to_str().unwrap())
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools view*******\n\n");
        if params.verbose {
            eprintln!("{}", String::from_utf8_lossy(&output.stderr));
            eprintln!("{}", "Samtools sort complete; Running samtools view");
        }
        info!("{}", String::from_utf8_lossy(&output.stderr));
        info!("{}", "Samtools sort complete; Running samtools view");
        let output = Command::new("samtools")
                        .arg("view")
                        .arg("-b")
                        .arg("-@")
                        .arg(params.threads.to_string())
                        .arg("-o")
                        .arg(align_sorted_bam.to_str().unwrap())
                        .arg(align_sorted_sam.to_str().unwrap())
                        .stderr(Stdio::piped())
                        .stdout(Stdio::piped())
                        .output()
                         .expect("\n\n*******Failed to execute samtools sort*******\n\n");
        if params.verbose {
            eprintln!("{}", String::from_utf8_lossy(&output.stderr));
            eprintln!("{}", "Samtools view complete; Running samtools index");
        }
        info!("{}", String::from_utf8_lossy(&output.stderr));
        info!("{}", "Samtools view complete; Running samtools index");
    }
  //   /app/software/CellRanger/6.0.1/lib/bin/STAR --genomeDir $transcriptome/star --readFilesIn ${fq3} --readNameSeparator space \
  // --runThreadN 24 --outSAMunmapped Within KeepPairs --outSAMtype BAM SortedByCoordinate

    if params.aligner == "STAR" {
        if params.verbose {
        eprintln!("{}", "Aligning reads using STAR");
        }
        info!("{}", "Aligning reads using STAR");
        let output = Command::new(&params.aligner_loc)
                        .arg("--outFileNamePrefix")
                        .arg(outfolder.to_str().unwrap())
                        .arg("--genomeDir")
                        .arg(params.genome.to_string())
                        .arg("--readFilesIn")
                        .arg(fastq.to_str().unwrap())
                        .arg("--readNameSeparator")
                        .arg("space")
                        .arg("--runThreadN")
                        .arg(params.threads.to_string())
                        .arg("--outSAMunmapped") 
                        .arg("Within")
                        .arg("KeepPairs") 
                        .arg("--outSAMtype") 
                        .arg("BAM")
                        .arg("SortedByCoordinate")
                        .arg("--outSAMattributes") 
                        .arg("All")
                        .arg("--readFilesCommand")
                        .arg("zcat")
                        .stderr(Stdio::piped())
                        .stdout(Stdio::piped())
                         .output()
                         .expect("\n\n*******Failed to execute STAR*******\n\n");
        if params.verbose {
            eprintln!("{}", String::from_utf8_lossy(&output.stderr));
            eprintln!("{}", "STAR complete; running samtools index");
        }
    }
    let output = Command::new("samtools")
                    .arg("index")
                    .arg("-@")
                    .arg(params.threads.to_string())
                    .arg(align_sorted_bam.to_str().unwrap())
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools index*******\n\n");
    if params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }
    info!("{}", String::from_utf8_lossy(&output.stderr));
    if !params.keep {
        fs::remove_file(align_sorted_sam.to_str().unwrap())?;
        fs::remove_file(align_sam.to_str().unwrap())?;
        fs::remove_file(fastq.to_str().unwrap())?;
    }
    Ok(())
}




pub fn writer_fn (count_vec: Vec<Vec<Vec<u8>>>, params: &Paramsm) -> Result<(), Box<dyn Error>> {
        let counts_path = &params.output_path.join("counts.txt.gz");
        let counts_file = counts_path.to_str().unwrap();
        info!("\n\n\tWriting counts to : '{}'\n", counts_file);
        if params.verbose{
            eprintln!("Writing counts to : '{}'\n", counts_file);
        }
        let f = File::create(counts_file)?;
        let mut gz = GzBuilder::new()
                        .filename(counts_file)
                        .write(f, Compression::default());
        for result in count_vec {
                for subresult in result {
                    gz.write_all(&subresult)?;
                }
        }
        gz.finish()?;
        Ok(())
}




fn remove_whitespace(s: &mut String) {
    s.retain(|c| !c.is_whitespace());
}



fn fastq(params: &Paramsm) -> Result<(), Box<dyn Error>>{
    let outfastq_temp = &params.output_path.join("mutcaller_R1.fq.gz").clone().to_owned();
    let outfastq = outfastq_temp.to_str().unwrap();
    let split = "|BARCODE=".to_string();
    let mut cbvec = lines_from_file(&params.bcs);
    cbvec.sort_unstable();
    let _zip = true;
    let mut total_count: usize = 0;
    let mut nfound_count: usize = 0;
    let mut mmcb_count: usize = 0;
    let split_at = &params.umi_len + &params.cb_len;
    // let sep: Vec::<u8> = params.name_sep.as_bytes().to_vec();

    let fastq1 = &params.fastq1;
    let fastq2 = &params.fastq2;
    let _counts = (0u64, 0u64);
    let path = Path::new(&outfastq);
    let _file = match File::create(&path) {
        Err(_why) => panic!("couldn't open {}", path.display()),
        Ok(file) => file,
    };
    // let mut writer = io::stdout();
    let f = File::create(&outfastq)?;
    let mut writer = GzBuilder::new()
                            .filename(outfastq)
                            .write(f, Compression::default());
    parse_path(Some(fastq1), |parser1| {
        parse_path(Some(fastq2), |parser2| {
            each_zipped(parser1, parser2, |rec1, rec2| {
                if rec1.is_some() & rec2.is_some(){
                    let r1 = &rec1.unwrap();
                    let r2 = &rec2.unwrap();
                    if r1.seq().contains(&b"N"[0]) | r2.seq().contains(&b"N"[0]){
                        nfound_count += 1;
                        total_count +=1;
                    }else{
                        total_count +=1;
                        let (barcode, _seq) = &r1.seq().split_at(split_at.into());
                        let (cb, _seq) = barcode.split_at(params.cb_len as usize);
                        match cbvec.binary_search(&std::str::from_utf8(cb).unwrap().to_string()) {
                            Ok(_u) => {
                                let mut readout = RefRecord::to_owned_record(&r2);
                                let _some_x = vec![b" "];
                                let mut new_header = std::str::from_utf8(&readout.head()).unwrap().to_string();
                                remove_whitespace(&mut new_header);
                                let _ = new_header.push_str(&split);
                                let _ = new_header.push_str(&std::str::from_utf8(&barcode).unwrap().to_string());
                                readout.head = new_header.as_bytes().to_vec();

                                let _ = readout.write(&mut writer);
                            }
                            Err(_e) => {
                                mmcb_count +=1;
                            }
                        }
                    }
                }
                (true, true)
            })
            .expect("Invalid record.");
        })
        .expect("Unknown format for file 2.");
    })
    .expect("Unknown format for file 1.");
    info!("\n\n\tTotal number of reads processed: {}, {} of these had Ns, {} of these had BC not in whitelist\n", total_count, nfound_count, mmcb_count);
    if params.verbose {
        eprintln!("Total number of reads processed: {}, {} of these had Ns, {} of these had BC not in whitelist\n", total_count, nfound_count, mmcb_count);
    }
    Ok(())
}



fn process_variant(ref_id: u32, start: u32)->bam::Region{
    let region = bam::Region::new(ref_id,start - 1,start - 1);
    return region;
}

fn count_variants_helper(params: &Paramsm, variant: Variant)-> Option<Vec<Vec<u8>>> {
        if variant.class.clone().unwrap() == VariantClass::SNV {
            return Some(count_variants_snv(params, variant))
        }
        warn!("\n\n\tVariant type {:?} not currently supported", variant.class.as_ref().unwrap());
       if params.verbose {
            eprintln!("\n\n\tVariant type {:?} not currently supported", variant.class.unwrap());
        }
        return None
}


#[allow(unused_comparisons)]
fn count_variants_snv(params: &Paramsm, variant: Variant) -> Vec<Vec<u8>>{
    // eprintln!("Processing using cb and umi in header");
    let split = "|BARCODE=".to_string();
    let ibam_temp = &params.output_path.join("Aligned.sortedByCoord.out.bam").clone().to_owned();
    let ibam = ibam_temp.to_str().unwrap();
    let mut total: usize = 0;
    let mut err: usize = 0;
    let seqname = variant.seq;
    let start = variant.start.parse::<u32>().unwrap();
    let vname = variant.name;
    let mut reader = bam::IndexedReader::build()
        .additional_threads(*&params.threads as u16)
        .from_path(ibam).unwrap();
    let mut seqnames = Vec::new();
    let mut _result = "";
    let query_nt = variant.query_nt.chars().nth(0);
    let header = reader.header().clone();
    let hdata = header.reference_names();
    for seq in hdata {
        seqnames.push(seq)
    }
    let mut data = Vec::new();
    let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
    let region = process_variant(ref_id as u32, start);
    for record in reader.fetch(&&region).unwrap(){
    // for record in reader.fetch_by(&&region, |record| record.mapq() >= 4 && (record.flag().all_bits(0 as u16) || record.flag().all_bits(16 as u16))).unwrap(){
        total+=1;
        let readheader = match str::from_utf8(record.as_ref().unwrap().name()) {
            Ok(v) => v,
            Err(e) => panic!("\n\n*******Invalid UTF-8 sequence: {}*******\n\n", e),
        };
        let cbumi = match readheader.split(&split).nth(1){
            Some(v) => v.to_string(),
            None => {
                err+=1;
                continue
            },
        };
        if cbumi.len() <= params.cb_len+1 {
            err+=1;
            continue
        }
        let (cb, umi) = cbumi.split_at((params.cb_len+1).into());
        for entry in record.as_ref().unwrap().alignment_entries().unwrap() {
            if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
                if region.start() == ref_pos {
                    if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
                        if ref_nt as char == record_nt as char {
                            _result = "ref";
                        } else if record_nt as char == query_nt.unwrap(){
                            _result = "query";
                        } else {
                            _result = "other";
                        }
                            data.push(format!("{} {} {} {} {} {}", &cb, &umi, seqname, ref_pos, vname, _result))
                        }
                    } else {
                        continue
                    }
            } else {
                continue
            }        }
    }
    info!("\n\n\tFound {} reads spanning this variant!\n\tNumbers of errors: {}", total, err);
    if params.verbose{
        eprintln!("Found {} reads spanning this variant!\n\tNumbers of errors: {}", total, err);
    }
    data.sort();
    let mut out_vec = Vec::new();
    let cdata = data.into_iter().dedup_with_count();
    for (count, record) in cdata {
       let count_str = record+&" ".to_owned()+&(count.to_string()+&"\n".to_owned());
        out_vec.push(count_str.as_bytes().to_owned());
    }
    return out_vec;
}


// fn count_variants_mm(params: &Paramsm, variant: Variant) -> Vec<Vec<u8>>{
//     eprintln!("Processing using cb and umi in BAM tags");
//     // let split = "|BARCODE=".to_string();
//     let ibam = "Aligned.mm2.bam";
//     let mut total: usize = 0;
//     let seqname = variant.seq;
//     let start = variant.start.parse::<u32>().unwrap();
//     let vname = variant.name;
//     let mut reader = bam::IndexedReader::build()
//         .additional_threads(*&params.threads as u16)
//         .from_path(ibam).unwrap();
//     let mut seqnames = Vec::new();
//     let mut cb;
//     let mut umi;
//     let mut result = "null";
//     let query_nt = variant.query_nt as char;
//     let header = reader.header().clone();
//     let hdata = header.reference_names();
//     for seq in hdata {
//         seqnames.push(seq)
//     }
//     let mut data = Vec::new();
//     let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
//     let region = process_variant(ref_id as u32, start);
//     for record in reader.fetch_by(&&region, |record| record.mapq() > 4 && (record.flag().all_bits(0 as u16) || record.flag().all_bits(16 as u16))).unwrap(){
//         total+=1;
//         match record.as_ref().unwrap().tags().get(b"CB") {
//             Some( bam::record::tags::TagValue::String(cba, _)) => {
//                 cb = str::from_utf8(&cba).unwrap().to_string();
//             },
//             _ => panic!("Unexpected type"),
//         }
//         match record.as_ref().unwrap().tags().get(b"UB") {
//             Some( bam::record::tags::TagValue::String(uba, _)) => {
//                 umi = str::from_utf8(&uba).unwrap().to_string();
//             },
//             _ => panic!("Unexpected type"),
//         }
//         for entry in record.as_ref().unwrap().alignment_entries().unwrap() {
//             if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
//                 if region.start() == ref_pos {

//                     if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
//                         if ref_nt as char == record_nt as char {
//                             result = "ref";
//                         } else if record_nt as char == query_nt{
//                             result = "query";
//                         } else {
//                             result = "other";
//                         }
//                             data.push(format!("{} {} {} {} {} {}", &cb, &umi, seqname, ref_pos, vname, result))
//                         }
//                     } else {
//                         continue
//                     }
//             } else {
//                 continue
//             }        }
//     }
//     eprintln!("Found {} reads spanning this variant!", total);
//     data.sort();
//     let mut out_vec = Vec::new();
//     let cdata = data.into_iter().dedup_with_count();
//     for (count, record) in cdata {
//         let count_str = record+&" ".to_owned()+&(count.to_string()+&"\n".to_owned());
//         // let count_str = record+&" ".to_owned()+&(count.to_string())+&"\n".to_owned();
//         out_vec.push(count_str.as_bytes().to_owned());
//     }
//     return out_vec;
// }


// fn count_star(params: &Paramsm) {
//     let ibam = "Aligned.mm2.bam";
//     let split = "|BARCODE=".to_string();
//     let joiner = "_".to_string();
//     eprintln!("Counting star reads");
//     let mut total: usize = 0;
//     let mut goodreadcount: usize = 0;
//     let (_read_threads, _write_threads) = if (*&params.threads as i8) > 2{
//         (((*&params.threads/2) -1) as usize, ((*&params.threads/2) -1) as usize)
//     } else {
//         (0 as usize, 0 as usize)
//     };

//     let reader = bam::BamReader::from_path(ibam.to_string(), 0).unwrap();
//     let _output = std::io::BufWriter::new(io::stdout());
//     let header = reader.header().clone();
//     let data = header.reference_names();
//     let mut seqnames = Vec::new();
//     for seq in data {
//         seqnames.push(seq)
//     }
//     for record in reader {
//         total += 1;
//         let newrecord = record.unwrap();
//         let seqname = match str::from_utf8(&newrecord.name()) {
//             Ok(v) => v,
//             Err(e) => panic!("\n\n*******Invalid UTF-8 sequence: {}*******\n\n", e),
//         };
//         let cbumi= seqname.split(&split).nth(1).unwrap().to_string();
//         let _modified_name = seqname.replace(&split, &joiner);
//         let (cb_umi_s1, cb_umi_s2) = cbumi.split_at((params.cb_len+1).into());
//         let mut good_read = false;
//         let cigarmatch = format!("{}M", *&params.read_len);
//         let cigar = newrecord.cigar().to_string();
//         if cigar == cigarmatch{
//             good_read = true
//         }
//         if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
//             goodreadcount += 1;
//             println!("{} {} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string(), newrecord.start());
//         }
//     }
//     eprintln!("Completed; {} total reads processed!", &total);
//     eprintln!("{} good reads counted!", &goodreadcount);
// }



// fn count_kallisto(params: &Paramsm) {
//     eprintln!("Counting kallisto reads");
//     let mut total: usize = 0;
//     let ibam = "Aligned.mm2.bam";
//     let split = "|BARCODE=".to_string();
//     let joiner = "_".to_string();
//     let mut goodreadcount: usize = 0;
//     let (_read_threads, _write_threads) = if (*&params.threads as i8) > 2{
//         (((*&params.threads/2) -1) as usize, ((*&params.threads/2) -1) as usize)
//     } else {
//         (0 as usize, 0 as usize)
//     };
//     let reader = bam::BamReader::from_path(ibam.to_string(), 0).unwrap();
//     let _output = io::BufWriter::new(io::stdout());
//     let header = reader.header().clone();
//     let data = header.reference_names();
//     let mut seqnames = Vec::new();
//     for seq in data {
//         seqnames.push(seq)
//     }
//     for record in reader {
//         total += 1;
//         let newrecord = record.unwrap();
//         let seqname = match str::from_utf8(&newrecord.name()) {
//             Ok(v) => v,
//             Err(e) => panic!("\n\n*******Invalid UTF-8 sequence: {}*******\n\n", e),
//         };
//         let cbumi= seqname.split(&split).nth(1).unwrap().to_string();
//         let _modified_name = seqname.replace(&split, &joiner);
//         let (cb_umi_s1, cb_umi_s2) = cbumi.split_at((params.cb_len+1).into());
//         let mut good_read = false;
//         let cigarmatch = format!("{}M", *&params.read_len);
//         let cigar = newrecord.cigar().to_string();
//         if cigar == cigarmatch{
//             good_read = true
//         }
//         if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
//             goodreadcount += 1;
//             println!("{} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string());
//         }
//     }
//     eprintln!("Completed; {} total alignments processed!", &total);
//     eprintln!("{} good alignments counted!", &goodreadcount);
// }

fn lines_from_file(filename: &str) -> Vec<String> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(_why) => panic!("\n\n*******couldn't open {}*******\n\n", path.display()),
        Ok(file) => file,
    };
    if path.extension() == Some(OsStr::new("gz")){
        let buf = BufReader::new(read::GzDecoder::new(file));
        buf.lines()
            .map(|l| l.expect("\n\n*******Could not parse line*******\n\n"))
            .collect()
    }else{
        let buf = BufReader::new(file);
        buf.lines()
            .map(|l| l.expect("\n\n*******Could not parse line*******\n\n"))
            .collect()
    }
}


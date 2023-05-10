/**


Full pipeline in 1 command...

cd ~/develop/mutCaller/tests
cargo build --release
../target/release/mutcaller --help

bc=~/develop/mutCaller/data/737K-august-2016.txt.gz

fa=/Users/sfurlan/refs/genome.fa
fa=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
../target/release/mutcaller \
                        -t 8 -g $fa -b $bc -v variants.tsv \
                        --fastq1 sequencer_R1.fastq.gz \
                        --fastq2 sequencer_R2.fastq.gz


*/


extern crate csv;
extern crate clap;
extern crate bam;
extern crate serde;
extern crate fastq;
extern crate itertools;


use std::io;
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
use simple_log::LogConfigBuilder;
use fastq::RefRecord;
use crate::mutcaller::fastq::Record;
// use crate::fastq::Record;
// use std::ffi::OsStr;
use flate2::{read};
use std::process::{Command, Stdio};
use std::ffi::OsStr;
use std::fs;


#[derive(Deserialize)]
struct Variant {
    seq: String,
    start: String,
    ref_nt: char,
    query_nt: char,
    name: String,
}

// Implement `Display` for `Variant`.
impl fmt::Display for Variant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use `self.number` to refer to each positional data point.
        write!(f, "seq: {} start: {} ref_nt: {} query_nt: {} name: {}", self.seq, self.start, self.ref_nt, self.query_nt, self.name)
    }
}


#[derive(Debug)]
struct Params {
    fastq1: String,
    fastq2: String,
    genome: String,
    bcs: String,
    umi_len: usize,
    cb_len: usize,
    threads: usize,
    aligner: String,
    variants: String,
    read_len: usize,
    output: String,
    keep: bool,
    verbose: bool,
}

// #[derive(Debug)]
// struct CountBamParams {
//     bam:
//     method:
//     stringsep:
//     cbsep:
//     umisep:
// }




fn load_params() -> Params {
    let yaml = load_yaml!("../cli.yml");
    let params = App::from_yaml(yaml).get_matches();
        let fastq1 = params.value_of("fastq1").unwrap();
        let fastq2 = params.value_of("fastq2").unwrap();
        let output = params.value_of("output").unwrap_or("counts_mm.txt.gz");
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
        let aligner = params.value_of("aligner").unwrap_or("mm2");
        let variantstring = params.value_of("variants").unwrap();
        let mut verbose = true;
        if params.is_present("quiet") {
                verbose = false
        };
        let mut keep = false;
        if params.is_present("keep_files") {
                keep = true
        };
        // let verbose = match params.value_of("verbose") {
        //     Some(val) => println!("true"),
        //     None => println!("false"),
        // };
        // let keep = match params.value_of("keep_files"){
        //     Some(val) => println!("true"),
        //     None => println!("false"),
        // };
        // let verbose = match params.value_of("verbose"){
        //     Err(..) => panic!("Verbose flag not recognized"),
        //     Ok(verbose) => verbose,
        // };
        // let keep = match params.value_of("keep_files"){
        //     Err(..) => panic!("keep_files flag not recognized"),
        //     Ok(verbose) => verbose,
        // };

        Params{
            fastq1: fastq1.to_string(),
            fastq2: fastq2.to_string(),
            genome: genome.to_string(),
            output: output.to_string(),
            bcs: bcs.to_string(),
            threads: threads as usize,
            umi_len: umi_len as usize,
            cb_len: cb_len as usize,
            aligner: aligner.to_string(),
            variants: variantstring.to_string(),
            read_len: read_len as usize,
            keep: keep,
            verbose: verbose,
        }
}




fn read_csv(params: &Params) -> Result<Vec<Variant>, Box<dyn Error>> {
    // Build the CSV reader and iterate over each record.
//     let data = "\
// seq\tstart\tref_nt\tquery_nt\tname
// chr12\t112450407\tA\tG\tPTPN11_227A>G
// chr2\t208248389\tG\tA\tIDH1_132G>A
// chr17\t7674220\tC\tT\tTP53_248C>T";
    eprintln!("Opening variants file: {}\n", &params.variants.to_string());
    let file = File::open(&params.variants.to_string()).unwrap();
    let reader = BufReader::new(file);
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(reader);
    let mut csvdata = Vec::new();
    for result in rdr.deserialize() {
        let record: Variant = result?;
        csvdata.push(record);
    }
    // let dummyvariant: Variant = Variant {seq: String::from("chr12"),
    //         start: String::from("112450407"),
    //         ref_nt: 'A',
    //         query_nt: 'G',
    //         name: String::from("PTPN11_227A>G"),
    //     };
    // csvdata.push(dummyvariant);
    Ok(csvdata)
}



pub fn mutcaller_run() {
    let config = LogConfigBuilder::builder()
        .path("./mutcaller.log")
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("debug")
        .output_file()
        .build();
    let _ = simple_log::new(config);
    info!("starting!");



    let params = load_params();
        if params.verbose {
        eprintln!("\n\n\n\nParsing Parameters!\n");
    }
    if params.verbose {
        eprintln!("\n\n\n\nChecking programs and parsing variants!\n");
    }
    let _prog_test_res = test_progs();
    let csvdata = read_csv(&params).unwrap();
    
    if params.verbose {
        for variant in &csvdata {
                eprintln!("\nCorrectly processed variant: {}", variant);
        }
    }
    if params.verbose {
        eprintln!("\n\n\n\nRunning with {} thread(s)!\n", &params.threads);
        // eprintln!("Params: {:?} ", &params);
        eprintln!("Processing fastqs:\n\t{}\n\t{}\n", &params.fastq1, &params.fastq2);
    }
    let _fqr = fastq(&params);
    info!("done!");
    let _ar = align(&params);
    if params.aligner == "mm2" {
        // let csvdata = read_csv(&params).unwrap();
        let mut count_vec = Vec::new();
        for variant in csvdata {
            eprintln!("\nProcessing variant: {}", variant);
            eprintln!("\nOpening bam: {}", &"Aligned.mm2.bam".to_string());
            count_vec.push(count_variants_mm2(&params, variant));
            
        }
        let _none = writer_fn(count_vec, params.output.to_string());
        eprintln!("\n\nDone!!");
        return;
    }
    if params.aligner == "kallisto"{
        count_kallisto(&params);
        return;
    }
    if params.aligner == "STAR"{
        count_star(&params);
        return;
    }
}

// minimap2 --MD -a $fa -t 8 mutcaller_R1.fq.gz -o Aligned.mm2.sam
// samtools sort -@ 8 -o Aligned.mm2.sorted.sam Aligned.mm2.sam
// samtools view -b -@ 8 -o Aligned.mm2.sorted.sam Aligned.mm2.sorted.bam
// samtools index -@ 8 Aligned.mm2.sorted.bam

fn test_progs () -> Result<(), Box<dyn Error>>{
    let _output = Command::new("minimap2")
                    .arg("-h")
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                     .output()
                     .expect("\n\n*******Failed to execute minimap2*******\n\n");
    let _output = Command::new("samtools")
                    .arg("-h")
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                     .output()
                     .expect("\n\n*******Failed to execute samtools*******\n\n");
    // eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    Ok(())
}


fn align (params: &Params)-> Result<(), Box<dyn Error>> {
    // let mm_cmd = "/Users/sfurlan/.local/bin/minimap2";
    // let mm_args = "-h";
    // let mm_args_pre = format!("-a {} -t {} mutcaller_R1.fq.gz | samtools sort -@ {} | samtools view -@ {} -o Aligned.mm2.bam", params.genome.to_string(), params.threads.to_string(), params.threads.to_string(), params.threads.to_string());
    // let mm_args_vec = mm_args_pre.split(" ").collect::<Vec<&str>>();
    // let mm_args = OsString::new();
    // let mm_args = OsString::from(mm_args_pre);
    // let mm_cmd = "/Users/sfurlan/.local/bin/minimap2";
    // let mm_cmd = "ls";
    // let st_cmd = format!("sammtools index -@ {} Aligned.mm2.bam", params.threads.to_string());
    // let output = Command::new("echo")
    //                  .arg("Hello world")
    //                  .output()
    //                  .expect("Failed to execute command");

    // eprintln!("{:?}", &mm_args);
    // let output = Command::new(mm_cmd)
    //                 .arg(mm_args)
    //                  .output()
    //                  .expect("Failed to execute minimap2");
    // let output = Command::new("/Users/sfurlan/.local/bin/minimap2")
    //                 .arg("-h")
    //                  .output()
    //                  .expect("Failed to execute minimap2");
    eprintln!("{}", "Aligning reads using minimap2");
    let output = Command::new("minimap2")
                    .arg("--MD")
                    .arg("-a")
                    .arg(params.genome.to_string())
                    .arg("-t")
                    .arg(params.threads.to_string())
                    .arg("mutcaller_R1.fq.gz")
                    .arg("-o")
                    .arg("Aligned.mm2.sam")
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                     .output()
                     .expect("\n\n*******Failed to execute minimap2*******\n\n");
    eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    eprintln!("{}", "Minimap2 complete; Running samtools sort");
    let output = Command::new("samtools")
                    .arg("sort")
                    .arg("-@")
                    .arg(params.threads.to_string())
                    .arg("-o")
                    .arg("Aligned.mm2.sorted.sam")
                    .arg("Aligned.mm2.sam")
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools view*******\n\n");
    eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    eprintln!("{}", "Samtools sort complete; Running samtools view");
    let output = Command::new("samtools")
                    .arg("view")
                    .arg("-b")
                    .arg("-@")
                    .arg(params.threads.to_string())
                    .arg("-o")
                    .arg("Aligned.mm2.sorted.bam")
                    .arg("Aligned.mm2.sorted.sam")
                    .stderr(Stdio::piped())
                    .stdout(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools sort*******\n\n");
    eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    eprintln!("{}", "Samtools view complete; Running samtools index");
    let output = Command::new("samtools")
                    .arg("index")
                    .arg("-@")
                    .arg(params.threads.to_string())
                    .arg("Aligned.mm2.sorted.bam")
                    .stdout(Stdio::piped())
                    .stderr(Stdio::piped())
                    .output()
                     .expect("\n\n*******Failed to execute samtools index*******\n\n");
    eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    if !params.keep {
        fs::remove_file("Aligned.mm2.sorted.sam")?;
        fs::remove_file("Aligned.mm2.sam")?;
        fs::remove_file("mutcaller_R1.fq.gz")?;
    }
    Ok(())
}



fn writer_fn (count_vec: Vec<Vec<Vec<u8>>>, fname: String) -> Result<(), Box<dyn Error>> {
        let f = File::create(fname)?;
        let mut gz = GzBuilder::new()
                        .filename("counts_mm.txt.gz")
                        .write(f, Compression::default());
        for result in count_vec {
            for line in result {
                gz.write_all(&line)?;
            }
        }
        gz.finish()?;
        Ok(())
}




fn remove_whitespace(s: &mut String) {
    s.retain(|c| !c.is_whitespace());
}



fn fastq(params: &Params) -> Result<(), Box<dyn Error>>{
    let split = "|BARCODE=".to_string();
    let outfastq = "mutcaller_R1.fq.gz".to_string();
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
    eprintln!("Total number of reads processed: {}, {} of these had Ns, {} of these had BC not in whitelist\n", total_count, nfound_count, mmcb_count);
    Ok(())
}



fn process_variant(ref_id: u32, start: u32)->bam::Region{
    let region = bam::Region::new(ref_id,start - 1,start - 1);
    return region;
}


fn count_variants_mm2(params: &Params, variant: Variant) -> Vec<Vec<u8>>{
    eprintln!("Processing using cb and umi in header");
    let split = "|BARCODE=".to_string();
    let ibam = "Aligned.mm2.sorted.bam";
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
    let query_nt = variant.query_nt as char;
    let header = reader.header().clone();
    let hdata = header.reference_names();
    for seq in hdata {
        seqnames.push(seq)
    }
    let mut data = Vec::new();
    let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
    let region = process_variant(ref_id as u32, start);
    for record in reader.fetch_by(&&region, |record| record.mapq() >= 4 && (record.flag().all_bits(0 as u16) || record.flag().all_bits(16 as u16))).unwrap(){
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
                        } else if record_nt as char == query_nt{
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
    eprintln!("Found {} reads spanning this variant!\n\tNumbers of errors: {}", total, err);
    data.sort();
    let mut out_vec = Vec::new();
    let cdata = data.into_iter().dedup_with_count();
    for (count, record) in cdata {
       let count_str = record+&" ".to_owned()+&(count.to_string()+&"\n".to_owned());
        out_vec.push(count_str.as_bytes().to_owned());
    }
    return out_vec;
}


// fn count_variants_mm(params: &Params, variant: Variant) -> Vec<Vec<u8>>{
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


fn count_star(params: &Params) {
    let ibam = "Aligned.mm2.bam";
    let split = "|BARCODE=".to_string();
    let joiner = "_".to_string();
    eprintln!("Counting star reads");
    let mut total: usize = 0;
    let mut goodreadcount: usize = 0;
    let (_read_threads, _write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as usize, ((*&params.threads/2) -1) as usize)
    } else {
        (0 as usize, 0 as usize)
    };

    let reader = bam::BamReader::from_path(ibam.to_string(), 0).unwrap();
    let _output = std::io::BufWriter::new(io::stdout());
    let header = reader.header().clone();
    let data = header.reference_names();
    let mut seqnames = Vec::new();
    for seq in data {
        seqnames.push(seq)
    }
    for record in reader {
        total += 1;
        let newrecord = record.unwrap();
        let seqname = match str::from_utf8(&newrecord.name()) {
            Ok(v) => v,
            Err(e) => panic!("\n\n*******Invalid UTF-8 sequence: {}*******\n\n", e),
        };
        let cbumi= seqname.split(&split).nth(1).unwrap().to_string();
        let _modified_name = seqname.replace(&split, &joiner);
        let (cb_umi_s1, cb_umi_s2) = cbumi.split_at((params.cb_len+1).into());
        let mut good_read = false;
        let cigarmatch = format!("{}M", *&params.read_len);
        let cigar = newrecord.cigar().to_string();
        if cigar == cigarmatch{
            good_read = true
        }
        if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
            goodreadcount += 1;
            println!("{} {} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string(), newrecord.start());
        }
    }
    eprintln!("Completed; {} total reads processed!", &total);
    eprintln!("{} good reads counted!", &goodreadcount);
}



fn count_kallisto(params: &Params) {
    eprintln!("Counting kallisto reads");
    let mut total: usize = 0;
    let ibam = "Aligned.mm2.bam";
    let split = "|BARCODE=".to_string();
    let joiner = "_".to_string();
    let mut goodreadcount: usize = 0;
    let (_read_threads, _write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as usize, ((*&params.threads/2) -1) as usize)
    } else {
        (0 as usize, 0 as usize)
    };
    let reader = bam::BamReader::from_path(ibam.to_string(), 0).unwrap();
    let _output = io::BufWriter::new(io::stdout());
    let header = reader.header().clone();
    let data = header.reference_names();
    let mut seqnames = Vec::new();
    for seq in data {
        seqnames.push(seq)
    }
    for record in reader {
        total += 1;
        let newrecord = record.unwrap();
        let seqname = match str::from_utf8(&newrecord.name()) {
            Ok(v) => v,
            Err(e) => panic!("\n\n*******Invalid UTF-8 sequence: {}*******\n\n", e),
        };
        let cbumi= seqname.split(&split).nth(1).unwrap().to_string();
        let _modified_name = seqname.replace(&split, &joiner);
        let (cb_umi_s1, cb_umi_s2) = cbumi.split_at((params.cb_len+1).into());
        let mut good_read = false;
        let cigarmatch = format!("{}M", *&params.read_len);
        let cigar = newrecord.cigar().to_string();
        if cigar == cigarmatch{
            good_read = true
        }
        if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
            goodreadcount += 1;
            println!("{} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string());
        }
    }
    eprintln!("Completed; {} total alignments processed!", &total);
    eprintln!("{} good alignments counted!", &goodreadcount);
}

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


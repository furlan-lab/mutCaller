/**
this module handles the ALIGNED functions in mutcaller
*/

extern crate simplelog;
extern crate csv;
extern crate clap;
extern crate bam;
extern crate serde;
extern crate fastq;
extern crate itertools;
extern crate rayon;


use clap::{App, load_yaml};
use std::{env, fs, path::Path, path::PathBuf};
use std::error::Error;
use csv::ReaderBuilder;
use std::fs::File;
use std::io::{Write, BufReader};
use flate2::{GzBuilder, Compression};
use std::time::{Instant};
#[cfg(not(feature = "paris"))]
use log::*;
use simplelog::{Config, WriteLogger, CombinedLogger, LevelFilter};
use crate::mutcaller::{Variant, count_variants_helper, classify_variant};
use crate::vcf::{guess_vcf, guess_compression, read_vcf_compressed, read_vcf_uncompressed};
use rayon::prelude::*;


#[derive(Debug)]
pub struct Params {
    pub bam: String,
    pub threads: usize,
    pub variants: String,
    pub output_path: Box<Path>,
    pub verbose: bool,
    pub umi_tag: String,
    pub cb_tag: String,
    pub vcf: bool,
    pub qual: f64,
}


pub fn load_params() -> Params {
    let yaml = load_yaml!("../cli.yml");
    let matches = App::from_yaml(yaml).get_matches();
    let mut _verbose = true;
    let countbam_params = matches.subcommand_matches("ALIGNED").unwrap();
    let bam = countbam_params.value_of("bam").unwrap();
    let output = countbam_params.value_of("output").unwrap_or("out");
    let threads = countbam_params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap();
    let variantstring = countbam_params.value_of("variants").unwrap();
    if countbam_params.is_present("quiet") {
            _verbose = false
    };
    let outpath = Path::new(&output);
    let cb_tag = countbam_params.value_of("cb").unwrap_or("CB").to_string();
    let umi_tag = countbam_params.value_of("umi").unwrap_or("XM").to_string();
    let qual = countbam_params.value_of("qual").unwrap_or("95.0");
    let qual = qual.to_string().parse::<f64>().unwrap();
    let params = Params{
        bam: bam.to_string(),
        output_path: outpath.into() ,
        threads: threads as usize,
        variants: variantstring.to_string(),
        verbose: _verbose,
        cb_tag: cb_tag,
        umi_tag: umi_tag,
        vcf: false,
        qual: qual
    };
    return check_params(params).unwrap()

}


pub fn check_params(params: Params) -> Result<Params, Box<dyn Error>>{
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
    let is_vcf = guess_vcf(&params.variants);
    CombinedLogger::init(vec![
        #[cfg(not(feature = "termcolor"))]
        WriteLogger::new(
            LevelFilter::Info,
            Config::default(),
            File::create(log_file).unwrap(),
        ),
    ])
    .unwrap();
    Ok(Params{
            bam: params.bam,
            threads: params.threads,
            output_path: params.output_path,
            variants: params.variants,
            verbose: params.verbose,
            cb_tag: params.cb_tag,
            umi_tag: params.umi_tag,
            vcf: is_vcf.unwrap(),
            qual: params.qual
    })
}


fn read_csv(params: &Params) -> Result<Vec<Variant>, Box<dyn Error>> {
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
    Ok(csvdata)
}


pub fn cleanup(filename: &Path, warn: bool) -> std::io::Result<()> {
    if Path::new(filename).exists(){
        fs::remove_file(filename.to_str().unwrap())?;
        Ok(())
    }else {
        if warn{
            warn!("\n\n\tFile does not exist: '{:?}'\n", filename);
        }
        Ok(())
    }
}


pub fn countbam_run() {
    let start = Instant::now();
    let params = load_params();
    let csvdata = {
        if params.vcf {
            let is_compressed = guess_compression(&params.variants);
            if is_compressed.unwrap() {
                read_vcf_compressed(&params.variants, &params.qual, &params.verbose)
            } else {
                read_vcf_uncompressed(&params.variants, &params.qual, &params.verbose)
            }
        } else {
            Ok(read_csv(&params).unwrap())
        }
    };
    let mut classified_variants = Vec::new();
    for variant in csvdata.as_ref().unwrap() {
        let classified_variant = classify_variant(&variant);
        if params.verbose {
            eprintln!("Correctly parsed and classified variant: {}\n\n", classified_variant.as_ref().unwrap());
        }
        info!("\tCorrectly parsed and classified variant: {}\n\n", classified_variant.as_ref().unwrap());
        classified_variants.push(classified_variant.unwrap());
    }

    info!("\n\n\tRunning with {} thread(s)!\n", &params.threads);
    if params.verbose {
        eprintln!("\n\nRunning with {} thread(s)!\n", &params.threads);
    }

    if params.threads==1 {
        let count_vec = count_helper(&params, classified_variants);
        let _none = writer_fn(count_vec, &params);
    } else {
        let mut count_vec = Vec::new();
        let variants = classified_variants;
        let chunk_iter = variants.par_chunks(params.threads)
                                .map(|vec_variants_chunk| count_helper(&params, vec_variants_chunk.to_vec()));
        chunk_iter.collect_into_vec(&mut count_vec);
    }
    
    let duration = start.elapsed();
    info!("\n\n\tDone!!\n");
    info!("\n\n\tTime elapsed is: {:?}\n", duration);
    if params.verbose{
        eprintln!("\n\nDone!!");
        eprintln!("\n\nTime elapsed is: {:?}", duration);
    }
    return;

}



fn count_helper (params: &Params, csvdata: Vec<Variant>) -> Vec<Vec<Vec<u8>>>{
    let mut count_vec = Vec::new();
    for variant in csvdata {
            let classified_variant = classify_variant(&variant);
            info!("\n\n\tProcessing variant: {}\n", variant);
            info!("\n\n\tOpening bam: {}\n", &params.bam);
            if params.verbose{
                eprintln!("\nProcessing variant: {}", variant);
                eprintln!("\nOpening bam: {}", &params.bam);
            }
            count_vec.push(count_variants_helper(None, Some(&params), classified_variant.unwrap()).unwrap());
        }
    return count_vec
}

pub fn writer_fn (count_vec: Vec<Vec<Vec<u8>>>, params: &Params) -> Result<(), Box<dyn Error>> {
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

pub fn get_current_working_dir() -> std::io::Result<PathBuf> {
    env::current_dir()
}


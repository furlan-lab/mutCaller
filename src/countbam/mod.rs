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
// use std::str;
use std::error::Error;
use csv::ReaderBuilder;
use std::fs::File;
use std::io::{Write, BufReader};
// use itertools::Itertools;
use flate2::{GzBuilder, Compression};
// use rayon::prelude::*;
use std::time::{Instant};
// use std::sync::{Arc, Mutex};
#[cfg(not(feature = "paris"))]
use log::*;
use simplelog::{Config, WriteLogger, CombinedLogger, LevelFilter};
use crate::mutcaller::{Variant, count_variants_helper, classify_variant};
use crate::vcf::{guess_vcf, guess_compression, read_vcf_compressed, read_vcf_uncompressed};

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
    let mut params = load_params();
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
    // for variant in &csvdata {
    //     if params.verbose {
    //         eprintln!("\nCorrectly parsed variant: {}", variant);
    //     }
    //     info!("\n\n\tCorrectly parsed variant: {}\n", variant);
    // }
    for variant in csvdata.as_ref().unwrap() {
        if params.verbose {
            eprintln!("\nCorrectly parsed variant: {}", variant);
        }
        info!("\n\n\tCorrectly parsed variant: {}\n", variant);
        // count_vec.push(count_variants(&params, variant.clone()));
        // for svariant in variant{
        //     count_vec.push(count_variants(&params, svariant.clone()));
        // }
    }

    info!("\n\n\tRunning with {} thread(s)!\n", &params.threads);
    if params.verbose {
        eprintln!("\n\nRunning with {} thread(s)!\n", &params.threads);
    }

    params.threads = 1;
    if params.threads==1 {

        let count_vec = count_helper_single(&params, csvdata.unwrap());
        let _none = writer_fn(count_vec, &params);

    } else {
        eprintln!("{:?}", "not prallelized yet");
        // rayon::ThreadPoolBuilder::new().num_threads(params.threads).build_global().unwrap();
        // let count_vec = count_helper_multiple(&params, csvdata.unwrap());
        // let _none = writer_fn(count_vec, &params);
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


// fn count_helper_multiple (params: &Params, csvdata: Vec<Variant>) -> Vec<Vec<Vec<u8>>>{
//     let count_vecs = Vec::new();
//     csvdata.into_par_iter().for_each(|variant| {
//         let classified_variant = classify_variant(&variant);
//         let count_vec = count_variants_helper(None, Some(&params), classified_variant.unwrap());
//         count_vecs.lock().unwrap().push(count_vec.unwrap());
//         // eprintln!("{}", &variant);
//     });
//     let out = Arc::try_unwrap(count_vec).unwrap();
//     let out = out.into_inner().unwrap();
//     return out
// }



fn count_helper_single (params: &Params, csvdata: Vec<Variant>) -> Vec<Vec<Vec<u8>>>{
    let mut count_vec = Vec::new();
    for variant in csvdata {
            info!("\n\n\tProcessing variant: {}\n", variant);
            info!("\n\n\tOpening bam: {}\n", &params.bam);
            if params.verbose{
                eprintln!("\nProcessing variant: {}", variant);
                eprintln!("\nOpening bam: {}", &params.bam);
            }
            let classified_variant = classify_variant(&variant);
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


// fn string_pop(slice: &[u8]) -> &[u8; 2] {
//     slice.try_into().expect("slice with incorrect length")
// }


// fn count_variants_wrapper(params: &Params, variant: Variant) -> Vec<Vec<u8>>{
//     // let split = "|BARCODE=".to_string();
//     let ibam = &params.bam;
//     let mut total: usize = 0;
//     let mut err: usize = 0;
//     let seqname = variant.seq;
//     let start = variant.start.parse::<u32>().unwrap();
//     let vname = variant.name;
//     let query_nt = variant.query_nt.chars().nth(0).unwrap();
//     let mut reader = bam::IndexedReader::build()
//         // .additional_threads(*&params.threads as u16)
//         .from_path(ibam).unwrap();
//     let mut seqnames = Vec::new();
//     let mut _result = "";
//     let header = reader.header().clone();
//     let hdata = header.reference_names();
//     for seq in hdata {
//         seqnames.push(seq)
//     }
//     let mut _cb = "NULL".to_string();
//     let mut _umi = "NULL".to_string();
//     let mut data = Vec::new();
//     let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
//     let region = process_variant(ref_id as u32, start);
//     let cb_tag_b = string_pop(params.cb_tag.as_bytes());
//     let umi_tag_b = string_pop(params.umi_tag.as_bytes());
//     for record in reader.fetch(&&region).unwrap(){
//         total+=1;
//         match record.as_ref().unwrap().tags().get(cb_tag_b) {
//             Some( bam::record::tags::TagValue::String(cba, _)) => {
//                 _cb = str::from_utf8(&cba).unwrap().to_string();
//             },
//             _ => {
//                 // eprintln!("ERROR: 'CB' not found");
//                 err+=1;
//                 continue
//             }
//         }
//         match record.as_ref().unwrap().tags().get(umi_tag_b) {
//             Some( bam::record::tags::TagValue::String(uma, _)) => {
//                 _umi = str::from_utf8(&uma).unwrap().to_string();
//             },
//             _ => {
//                 // eprintln!("ERROR: 'CB' not found");
//                 err+=1;
//                 continue
//             }
//         }
//          for entry in record.as_ref().unwrap().alignment_entries().unwrap(){
//            if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
//                 if region.start() == ref_pos {
//                     if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
//                         if ref_nt as char == record_nt as char {
//                             _result = "ref";
//                         } else if record_nt as char == query_nt{
//                             _result = "query";
//                         } else {
//                             _result = "other";
//                         }
//                             // eprintln!("pushting data");
//                             data.push(format!("{} {} {} {} {} {}", _cb, _umi, seqname, ref_pos, vname, _result))
//                         }
//                     } else {
//                         continue
//                     }
//             } else {
//                 continue
//             }        
//         }
//     }

//     info!("\n\n\tFound {} reads spanning this variant!\n\tNumbers of errors: {}\n", total, err);
//     if params.verbose{
//         eprintln!("Found {} reads spanning this variant!\n\tNumbers of errors: {}\n", total, err);
//     }
    
//     data.sort();
//     let mut out_vec = Vec::new();
//     let cdata = data.into_iter().dedup_with_count();
//     for (count, record) in cdata {
//        let count_str = record+&" ".to_owned()+&(count.to_string()+&"\n".to_owned());
//         out_vec.push(count_str.as_bytes().to_owned());
//     }
//     return out_vec;
// }



pub fn get_current_working_dir() -> std::io::Result<PathBuf> {
    env::current_dir()
}


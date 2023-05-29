
/**


Full pipeline in 1 command...

cd ~/develop/mutCaller
cargo build --release
target/release/mutcaller --help

bc=~/develop/mutCaller/data/737K-august-2016.txt.gz

fa=/Users/sfurlan/refs/genome.fa
fa=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
../target/release/mutcaller \
                        -t 8 -g $fa -b $bc -v variants.tsv \
                        --fastq1 sequencer_R1.fastq.gz \
                        --fastq2 sequencer_R2.fastq.gz



cd /fh/scratch/delete90/furlan_s/targ_reseq/230117_Sami/AML_1101_merge/align_nodedup
ml SAMtools/1.11-GCC-10.2.0
samtools view -bs 42.01 aln_sncr_fc_Malig.mdtag.bam > subsampleMalig.mdtag.bam

cd /fh/scratch/delete90/furlan_s/targ_reseq/230117_Sami/AML_1101_merge/align_nodedup

~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v ~/develop/mutCaller/tests/variants.tsv -t 1 -o test1
~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v ~/develop/mutCaller/tests/variants.tsv -t 4


~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 1
~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 22

~/develop/mutCaller/target/release/mutcaller ALIGNED -t 20 -b aln_sncr_fc_Malig.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv



~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 1
~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 22

sbatch -n 1 -c 1 -p campus-new -M gizmo --mem-per-cpu=21000MB --wrap='~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 1 -o out1'
sbatch -n 1 -c 30 -p campus-new -M gizmo --mem-per-cpu=21000MB --wrap='~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 30 -o out30'
sbatch -n 1 -c 1 -p campus-new -M gizmo --mem-per-cpu=21000MB --wrap='~/develop/mutCaller/target/release/mutcaller ALIGNED -b aln_sncr_fc_Malig.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 1 -o out1_real'


head /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv


cat > breaksit.csv << EOL
seq,start,ref_nt,query_nt,name
chr14,31113320,T,C,HECTD1_25831T>C
EOL


sed -E 's/("([^"]*)")?,/\2\t/g' breaksit.csv > breaksit.tsv


~/develop/mutCaller/target/release/mutcaller ALIGNED -b aln_sncr_fc_Malig.mdtag.bam -v breaksit.tsv -t 1 -o out1_break

samtools view aln_sncr_fc_Malig.mdtag.bam | grep 'm64272e_230122_062448/100729383/ccs/9569_11265' > breaksit.bam
samtools view -b aln_sncr_fc_Malig.mdtag.bam "chr14:31113320-31113321" > breaksit2.bam


~/develop/mutCaller/target/release/mutcaller ALIGNED -b aln_sncr_fc_Malig.mdtag.bam -v ~/develop/mutCaller/tests/variants.tsv -t 1 -o test100
~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 1 -o out200


~/develop/mutCaller/target/release/mutcaller ALIGNED -b breaksit2.bam -v breaksit.tsv -t 1 -out_long


*/

extern crate simplelog;
extern crate csv;
extern crate clap;
extern crate bam;
extern crate serde;
extern crate fastq;
extern crate itertools;
extern crate rayon;


// use std::io;
use clap::{App, load_yaml};
use std::{env, str, fs, path::Path, path::PathBuf};
use std::error::Error;
use serde::Deserialize;
use std::fmt; 
use csv::ReaderBuilder;
use std::fs::File;
use std::io::{Write, BufReader};
use itertools::Itertools;
use flate2::GzBuilder;
use flate2::Compression;
use rayon::prelude::*;
use std::ops::ControlFlow;
// use simplelog::{Config, WriteLogger, SimpleLogger, CombinedLogger, LevelFilter, info, warn};
// use std::sync::Mutex;
// use core::iter::Map;
use std::time::{Instant};
// use std::path::Path;
// use fastq::parse_path;
// use fastq::each_zipped;
// use simple_log::LogConfigBuilder;
// use simple_log::{info, warn};
use std::sync::{Arc, Mutex};
// use bam::record::tags::TagValue;
// use crate::fastq::Record;
// use std::ffi::OsStr;
// use flate2::{read};
// use std::process::{Command, Stdio};
// use std::ffi::OsStr;
// use std::fs;
#[cfg(not(feature = "paris"))]
use log::*;
// use simplelog::*;
use simplelog::{Config, WriteLogger, CombinedLogger, LevelFilter};


#[derive(Deserialize)]
#[derive(Debug)]
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
pub struct Params {
    pub bam: String,
    pub threads: usize,
    pub variants: String,
    pub output_path: Box<Path>,
    pub verbose: bool,
    pub umi_tag: String,
    pub cb_tag: String,
}

// #[derive(Debug)]
// struct CountBamParams {
//     bam:
//     method:
//     stringsep:
//     cbsep:
//     umisep:
// }




pub fn load_params() -> Params {
    let yaml = load_yaml!("../cli.yml");
    let matches = App::from_yaml(yaml).get_matches();
    let mut _verbose = true;
    let countbam_params = matches.subcommand_matches("ALIGNED").unwrap();
    let bam = countbam_params.value_of("bam").unwrap();
    let output = countbam_params.value_of("output").unwrap_or("out");
    let threads = countbam_params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap();
    // let variantstring = countbam_params.value_of("variants").unwrap_or("/Users/sfurlan/develop/mutCaller/tests/variants.tsv");
    let variantstring = countbam_params.value_of("variants").unwrap();
    if countbam_params.is_present("quiet") {
            _verbose = false
    };
    let outpath = Path::new(&output);
    let cb_tag = countbam_params.value_of("cb").unwrap_or("CB").to_string();
    let umi_tag = countbam_params.value_of("umi").unwrap_or("XM").to_string();
    let params = Params{
        bam: bam.to_string(),
        output_path: outpath.into() ,
        threads: threads as usize,
        variants: variantstring.to_string(),
        verbose: _verbose,
        cb_tag: cb_tag,
        umi_tag: umi_tag,
    };
    return check_params(params).unwrap()

}


pub fn check_params(params: Params) -> Result<Params, Box<dyn Error>>{
    let _cu = cleanup(&params.output_path.join("mutcaller.log"), false);
    // let config = LogConfigBuilder::builder()
    //     .path(params.output_path.join("mutcaller.log").to_str().unwrap())
    //     .size(1 * 100)
    //     .roll_count(10)
    //     .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
    //     .level("info")
    //     .output_file()
    //     .build();
    // let _ = simple_log::new(config);
    let log_file_path = &params.output_path.join("mutcaller.log");
    let log_file = log_file_path.to_str().unwrap();
    // CombinedLogger::init(
    //     vec![
    //         TermLogger::new(LevelFilter::Warn, Config::default(), TerminalMode::Mixed, ColorChoice::Auto),
    //         WriteLogger::new(LevelFilter::Info, Config::default(), File::create(log_file).unwrap()),
    //     ]
    // ).unwrap();

    // CombinedLogger::init(vec![
    //     // #[cfg(feature = "termcolor")]
    //     // TermLogger::new(
    //     //     LevelFilter::Warn,
    //     //     Config::default(),
    //     //     TerminalMode::Mixed,
    //     //     ColorChoice::Always,
    //     // ),
    //     #[cfg(not(feature = "termcolor"))]
    //     // SimpleLogger::new(LevelFilter::Warn, Config::default()),
    //     WriteLogger::new(
    //         LevelFilter::Info,
    //         Config::default(),
    //         File::create(log_file).unwrap(),
    //     ),
    // ])
    // .unwrap();

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
        // #[cfg(feature = "termcolor")]
        // TermLogger::new(
        //     LevelFilter::Warn,
        //     Config::default(),
        //     TerminalMode::Mixed,
        //     ColorChoice::Always,
        // ),
        #[cfg(not(feature = "termcolor"))]
        // SimpleLogger::new(LevelFilter::Warn, Config::default()),
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
    })
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
    // let _cu = cleanup(&params.output_path.join("mutcaller.log"), false);
    // let config = LogConfigBuilder::builder()
    //     .path(&*params.output_path.join("mutcaller.log").to_str().unwrap())
    //     .size(1 * 100)
    //     .roll_count(10)
    //     .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
    //     .level("info")
    //     .output_file()
    //     .build();
    // let _ = simple_log::new(config);
    // info!("starting!");

    let csvdata = read_csv(&params).unwrap();
    

    for variant in &csvdata {
        if params.verbose {
            eprintln!("\nCorrectly parsed variant: {}", variant);
        }
        info!("\n\n\tCorrectly parsed variant: {}\n", variant);
    }

    info!("\n\n\tRunning with {} thread(s)!\n", &params.threads);
    if params.verbose {
        eprintln!("\n\nRunning with {} thread(s)!\n", &params.threads);
        // eprintln!("Params: {:?} ", &params);
    }


    // let mut count_vec: Vec<Vec<u8>> = Vec::new();
    // let mut count_vec = Vec::new();
    
    // ifparams.threads==1 {
    //     for variant in csvdata {
    //         eprintln!("\nProcessing variant: {}", variant);
    //         eprintln!("\nOpening bam: {}", &params.bam);
    //         count_vec.push(count_variants_mm2(&params, variant));
            
    //     }
    // }

    // for variant in csvdata {
    //     eprintln!("\nProcessing variant: {}", variant);
    //     eprintln!("\nOpening bam: {}", &params.bam);
    //     count_vec.push(count_variants_wrapper(&params, variant));
    // }

    // let count_vec = Mutex::new(count_vec);
    params.threads = 1;
    if params.threads==1 {

        let count_vec = count_helper_single(&params, csvdata);
        let _none = writer_fn(count_vec, &params);

    } else {
        rayon::ThreadPoolBuilder::new().num_threads(params.threads).build_global().unwrap();
        let count_vec = count_helper_multiple(&params, csvdata);
        let _none = writer_fn(count_vec, &params);
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


fn count_helper_multiple (params: &Params, csvdata: Vec<Variant>) -> Vec<Vec<Vec<u8>>>{
    let count_vecs = Arc::new(Mutex::new(Vec::new()));
    csvdata.into_par_iter().for_each(|variant| {
        let count_vec = count_variants_wrapper(&params, variant);
        count_vecs.lock().unwrap().push(count_vec);
        // eprintln!("{}", &variant);
    });
    let out = Arc::try_unwrap(count_vecs).unwrap();
    let out = out.into_inner().unwrap();
    return out
}



fn count_helper_single (params: &Params, csvdata: Vec<Variant>) -> Vec<Vec<Vec<u8>>>{
    let mut count_vec = Vec::new();
    for variant in csvdata {
            info!("\n\n\tProcessing variant: {}\n", variant);
            info!("\n\n\tOpening bam: {}\n", &params.bam);
            if params.verbose{
                eprintln!("\nProcessing variant: {}", variant);
                eprintln!("\nOpening bam: {}", &params.bam);
            }
            count_vec.push(count_variants_wrapper(&params, variant));
        }
    return count_vec
}

// fn writer_fn<T> (count_vec: Mutex<Vec<T>>, fname: String) -> Result<(), Box<dyn Error>> {
//         let f = File::create(&fname)?;
//         let mut gz = GzBuilder::new()
//                         .filename(fname)
//                         .write(f, Compression::default());
//         for result in count_vec {
//             for line in result {
//                 gz.write_all(&line)?;
//             }
//         }
//         gz.finish()?;
//         Ok(())
// }




// fn writer_fn (count_vec: Vec<Vec<Vec<u8>>>, fname: String) -> Result<(), Box<dyn Error>> {
//         let f = File::create(&fname)?;
//         let mut gz = GzBuilder::new()
//                         .filename(fname)
//                         .write(f, Compression::default());
//         for result in count_vec {
//             for line in result {
//                 gz.write_all(&line)?;
//             }
//         }
//         gz.finish()?;
//         Ok(())
// }

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



// fn remove_whitespace(s: &mut String) {
//     s.retain(|c| !c.is_whitespace());
// }



fn process_variant(ref_id: u32, start: u32)->bam::Region{
    let region = bam::Region::new(ref_id,start - 1,start - 1);
    return region;
}

fn string_pop(slice: &[u8]) -> &[u8; 2] {
    slice.try_into().expect("slice with incorrect length")
}

// fn test_md_tag(record: &Result<bam::record::Record, io::Error>) -> Result<(), io::Error> {
//     record.as_ref().unwrap().alignment_entries().unwrap();
//     Ok(())
// }

// fn test_md_tag (res: &Result<bam::record::Record, io::Error>) -> Result<(), io::Error> {
//     match res.as_ref().unwrap().alignment_entries(){
//         Ok(found) => {
//             // eprintln!("sequence ok")
//             eprintln!("{:?}", found);
//         }
//         Err(_e) => {eprintln!("sequence bad")}
//     }
//     Ok(())
// }

// fn find_error(id: &Id) -> Result<Item, String> {
//     Err(format!("Not found: {:?}", id))
// }

fn count_variants_wrapper(params: &Params, variant: Variant) -> Vec<Vec<u8>>{
    // let split = "|BARCODE=".to_string();
    let ibam = &params.bam;
    let mut total: usize = 0;
    let mut err: usize = 0;
    let seqname = variant.seq;
    let start = variant.start.parse::<u32>().unwrap();
    let vname = variant.name;
    let query_nt = variant.query_nt as char;
    let mut reader = bam::IndexedReader::build()
        // .additional_threads(*&params.threads as u16)
        .from_path(ibam).unwrap();
    let mut seqnames = Vec::new();
    let mut _result = "";
    let header = reader.header().clone();
    let hdata = header.reference_names();
    for seq in hdata {
        seqnames.push(seq)
    }
    // let mut cb = "NULL".to_string().as_bytes().to_vec();
    // let mut umi = "NULL".to_string().as_bytes().to_vec();
    let mut _cb = "NULL".to_string();
    let mut _umi = "NULL".to_string();
    let mut data = Vec::new();
    let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
    let region = process_variant(ref_id as u32, start);
    let cb_tag_b = string_pop(params.cb_tag.as_bytes());
    let umi_tag_b = string_pop(params.umi_tag.as_bytes());
    // for record in reader.fetch_by(&&region, |record| record.mapq() >= 0).unwrap(){
    for record in reader.fetch(&&region).unwrap(){
        total+=1;
        match record.as_ref().unwrap().tags().get(cb_tag_b) {
            Some( bam::record::tags::TagValue::String(cba, _)) => {
                _cb = str::from_utf8(&cba).unwrap().to_string();
            },
            _ => {
                // eprintln!("ERROR: 'CB' not found");
                err+=1;
                continue
            }
        }
        match record.as_ref().unwrap().tags().get(umi_tag_b) {
            Some( bam::record::tags::TagValue::String(uma, _)) => {
                _umi = str::from_utf8(&uma).unwrap().to_string();
            },
            _ => {
                // eprintln!("ERROR: 'CB' not found");
                err+=1;
                continue
            }
        }
        // eprintln!("{:?}", record);
        // ## try to unwrap
        // let do_steps = || -> Result<(), io::Error> {
        //     eprintln!("Testing read: {:?}", str::from_utf8(record.as_ref().unwrap().name()).unwrap());
        //     test_md_tag(&record)?;
        //     eprintln!("Alignment okay!");
        //     Ok(())
        // };

        // if let Err(_err) = do_steps() {
        //     if params.verbose{
        //         eprintln!("Failed to get alignment from read: {:?}", record.as_ref().unwrap().name());
        //     }
        //     error!("Failed to get alignment from read: {:?}", record.as_ref().unwrap().name());
        //     continue
        // }
        // eprintln!("Testing record: {:?}", str::from_utf8(record.as_ref().unwrap().name()).unwrap());

        // match record.as_ref().unwrap().tags().get(b"MD") {
        //     Some(bam::record::tags::TagValue::String(tag, _)) => {
        //         eprintln!("correct md tag" );
        //     },
        //     Some(_) => {
        //         eprintln!("incorrect md tag" );
        //         // &[0]
        //         },
        //     None => {
        //         eprintln!("No sequence found" );
        //         // &[0]
        //     },
        // };
        // eprintln!("{:?}", md_tag);



        // // maybe will work....
        // let _r = record.as_ref().unwrap().alignment_entries().unwrap().try_for_each(|entry| {
        //     // eprintln!("{:?}", entry.record_pos().unwrap());
        //     if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
        //      // eprintln!("{:?}", region.start());
        //         if region.start() == ref_pos {
        //             if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
        //                 if ref_nt as char == record_nt as char {
        //                     eprintln!("Ref found");
        //                     _result = "ref";
        //                 } else if record_nt as char == query_nt{
        //                     eprintln!("Query found");
        //                     _result = "query";
        //                 } else {
        //                     eprintln!("Other found");
        //                     _result = "other";
        //                 }   
        //                     eprintln!("pushting data");
        //                     data.push(format!("{} {} {} {} {} {}", _cb, _umi, seqname, ref_pos, vname, _result));
        //                     // return ControlFlow::Continue(())
        //                     return ControlFlow::Break(entry)
        //                 }
        //             return ControlFlow::Continue(())
        //             } else {
        //                 return ControlFlow::Continue(())
        //             }
        //         } else {
        //             return ControlFlow::Continue(())
        //         }        
        // });

        let iter = record.as_ref().unwrap().alignment_entries().unwrap();
        let mut new_iter = iter.map_while(|entry| {
            eprintln!("{:?}", entry);
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
                            eprintln!("pushting data");
                            data.push(format!("{} {} {} {} {} {}", _cb, _umi, seqname, ref_pos, vname, _result));
                            return Some(entry)
                    } else {
                        return Some(entry)
                    }

                } else {
                    return Some(entry)
                }
            } else {
                return Some(entry)
            }      
        });

        // for entry in record.as_ref().unwrap().alignment_entries(){
        //     // let do_steps = || -> Result<(), io::Error> {
        //     //     // eprintln!("Testing entry: {:?}", &entry);
        //     //     let rec = entry.record_pos_nt().unwrap();
        //     //     eprintln!("Rec pos: {:?}", char::from_u32(rec.0));
        //     //     eprintln!("Rec NT: {:?}", rec.1 as char);
        //     //     let refe = entry.ref_pos_nt().unwrap();
        //     //     eprintln!("Ref pos: {:?}", char::from_u32(refe.0));
        //     //     eprintln!("Ref NT: {:?}", refe.1 as char);
        //     //     // test_alignment(&record)?;
        //     //     // eprintln!("Entry okay!");
        //     //     Ok(())
        //     // };
        //     // if let Err(_err) = do_steps() {
        //     //     if params.verbose{
        //     //         eprintln!("Failed to get entry: {:?}", record.as_ref().unwrap().name());
        //     //     }
        //     //     error!("Failed to get entry: {:?}", record.as_ref().unwrap().name());
        //     //     continue
        //     // }
        //      // eprintln!("{:?}", entry);
        // //    if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
        // //         if region.start() == ref_pos {
        // //             if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
        // //                 if ref_nt as char == record_nt as char {
        // //                     _result = "ref";
        // //                 } else if record_nt as char == query_nt{
        // //                     _result = "query";
        // //                 } else {
        // //                     _result = "other";
        // //                 }
        // //                     eprintln!("pushting data");
        // //                     data.push(format!("{} {} {} {} {} {}", _cb, _umi, seqname, ref_pos, vname, _result))
        // //                 }
        // //             } else {
        // //                 continue
        // //             }
        // //     } else {
        // //         continue
        // //     }        
        // }

    }
     // data.push(format!("whatever"));
    info!("\n\n\tFound {} reads spanning this variant!\n\tNumbers of errors: {}\n", total, err);
    if params.verbose{
        eprintln!("Found {} reads spanning this variant!\n\tNumbers of errors: {}\n", total, err);
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



pub fn get_current_working_dir() -> std::io::Result<PathBuf> {
    env::current_dir()
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


// fn count_star(params: &Params) {
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



// fn count_kallisto(params: &Params) {
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

// fn lines_from_file(filename: &str) -> Vec<String> {
//     let path = Path::new(filename);
//     let file = match File::open(&path) {
//         Err(_why) => panic!("\n\n*******couldn't open {}*******\n\n", path.display()),
//         Ok(file) => file,
//     };
//     if path.extension() == Some(OsStr::new("gz")){
//         let buf = BufReader::new(read::GzDecoder::new(file));
//         buf.lines()
//             .map(|l| l.expect("\n\n*******Could not parse line*******\n\n"))
//             .collect()
//     }else{
//         let buf = BufReader::new(file);
//         buf.lines()
//             .map(|l| l.expect("\n\n*******Could not parse line*******\n\n"))
//             .collect()
//     }
// }


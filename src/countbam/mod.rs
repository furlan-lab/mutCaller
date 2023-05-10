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


// use std::io;
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
// use fastq::parse_path;
// use fastq::each_zipped;
use simple_log::LogConfigBuilder;
use bam::record::tags::TagValue;
// use crate::fastq::Record;
// use std::ffi::OsStr;
use flate2::{read};
// use std::process::{Command, Stdio};
use std::ffi::OsStr;
// use std::fs;


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
    bam: String,
    threads: usize,
    variants: String,
    output: String,
    verbose: bool,
    umi_tag: String,
    cb_tag: String,
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
        let bam = params.value_of("bam").unwrap();
        let output = params.value_of("output").unwrap_or("counts_mm.txt.gz");
        let threads = params.value_of("threads").unwrap_or("1");
        let threads = threads.to_string().parse::<usize>().unwrap();
        let variantstring = params.value_of("variants").unwrap();
        let mut verbose = true;
        if params.is_present("quiet") {
                verbose = false
        };
        let cb_tag = params.value_of("cb").unwrap_or("CB").to_string();
        let umi_tag = params.value_of("umi").unwrap_or("XM").to_string();

        Params{
            bam: bam.to_string(),
            output: output.to_string(),
            threads: threads as usize,
            variants: variantstring.to_string(),
            verbose: verbose,
            cb_tag: cb_tag,
            umi_tag: umi_tag,
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



pub fn countbam_run() {
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
        eprintln!("\n\n\n\nParsing variants!\n");
    }
    let csvdata = read_csv(&params).unwrap();
    
    if params.verbose {
        for variant in &csvdata {
                eprintln!("\nCorrectly processed variant: {}", variant);
        }
    }
    if params.verbose {
        eprintln!("\n\n\n\nRunning with {} thread(s)!\n", &params.threads);
        // eprintln!("Params: {:?} ", &params);
    }


    let mut count_vec = Vec::new();
    for variant in csvdata {
        eprintln!("\nProcessing variant: {}", variant);
        eprintln!("\nOpening bam: {}", &params.bam);
        count_vec.push(count_variants_mm2(&params, variant));
        
    }
    let _none = writer_fn(count_vec, params.output.to_string());
    eprintln!("\n\nDone!!");
    return;
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




// fn remove_whitespace(s: &mut String) {
//     s.retain(|c| !c.is_whitespace());
// }



fn process_variant(ref_id: u32, start: u32)->bam::Region{
    let region = bam::Region::new(ref_id,start - 1,start - 1);
    return region;
}


fn count_variants_mm2(params: &Params, variant: Variant) -> Vec<Vec<u8>>{
    eprintln!("Processing using cb and umi in header");
    // let split = "|BARCODE=".to_string();
    let ibam = &params.bam;
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
    // let mut cb = "NULL".to_string().as_bytes().to_vec();
    // let mut umi = "NULL".to_string().as_bytes().to_vec();
    let mut cb = "NULL".to_string();
    let mut umi = "NULL".to_string();
    let mut data = Vec::new();
    let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
    let region = process_variant(ref_id as u32, start);
    for record in reader.fetch_by(&&region, |record| record.mapq() >= 4 && (record.flag().all_bits(0 as u16) || record.flag().all_bits(16 as u16))).unwrap(){
        total+=1;
        match record.as_ref().unwrap().tags().get(b"CB") {
            // Some(TagValue::String(array_view, _)) => {
            //     cb = array_view.to_vec();
            // },
            Some(TagValue::Char(value)) => {
                cb= value.to_string();
            },
            _ => {
                // eprintln!("ERROR: 'CB' not found");
                err+=1;
                continue
            }
        }
        match record.as_ref().unwrap().tags().get(b"XM") {
            // Some(TagValue::String(array_view, _)) => {
            //     umi = array_view.to_vec();
            // },
            Some(TagValue::Char(value)) => {
                umi= value.to_string();
            },
            _ => {
                // eprintln!("ERROR: 'CB' not found");
                err+=1;
                continue
            }
        }
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
                            data.push(format!("{} {} {} {} {} {}", cb, umi, seqname, ref_pos, vname, _result))
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


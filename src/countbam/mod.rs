/**



cd ~/develop/mutCaller/tests
../target/release/countbam -b jaffa.ss.bam -w 0 -r
    
**/


extern crate csv;
extern crate clap;
extern crate bam;
extern crate serde;
extern crate itertools;


use clap::{App, load_yaml};
use std::str;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use flate2::GzBuilder;
use flate2::Compression;
use simple_log::info;
use simple_log::LogConfigBuilder;
// use crate::itertools::Itertools;




#[derive(Debug)]
struct Params {
    bam: String,
    threads: u16,
    method: String,
    stringsep: String,
    cbsep: String,
    umisep: String,
    sampleheader: i32,
    returnseq: bool,
    cbpos: usize,
    umipos: usize,
    outfile: String,
}

    // #[derive(Debug)]
    // struct Output {
    //     cb: String,
    //     umi: String,
    //     seq: String,
    //     start: i32,
    //     mapq: i32,
    //     seq: String,
    // }





fn load_params() -> Params {
    let yaml = load_yaml!("../cli.yml");
    let params = App::from_yaml(yaml).get_matches();
    let bam = params.value_of("bam").unwrap();
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap();
    let stringsep = params.value_of("stringsep").unwrap_or(";");
    let cbsep = params.value_of("cbsep").unwrap_or("XC=");
    let umisep = params.value_of("umisep").unwrap_or("XM=");
    let method = params.value_of("method").unwrap_or("header");
    let mut sampleheader = -1 as i32;
    if params.is_present("sampleheader") {
        sampleheader = params.value_of("sampleheader").unwrap().parse::<i32>().unwrap();
    };
    let mut returnseq = false;
    if params.is_present("returnseq") {
        returnseq = true;
    }
    let mut cbpos = 5 as usize;
    if params.is_present("cbpos") {
        cbpos = params.value_of("cbpos").unwrap().parse::<usize>().unwrap() - 1;
    }
    let mut umipos = 4 as usize;
    if params.is_present("umipos") {
        umipos = params.value_of("umipos").unwrap().parse::<usize>().unwrap() - 1;
    }
    if params.is_present("returnseq") {
        returnseq = true;
    }
    let outfile = params.value_of("outfile").unwrap_or("counts.txt.gz");
    Params {
        bam: bam.to_string(),
        threads: threads as u16,
        method: method.to_string(),
        stringsep: stringsep.to_string(),
        cbsep: cbsep.to_string(),
        umisep: umisep.to_string(),
        sampleheader: sampleheader as i32,
        returnseq: returnseq,
        cbpos: cbpos,
        umipos: umipos,
        outfile: outfile.to_string(),
    }
}




pub fn countbam_run() {
    let verbose = true;
    let config = LogConfigBuilder::builder()
        .path("./countbam.log")
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("debug")
        .output_file()
        .build();
    let _ = simple_log::new(config);
    info!("starting!");
    let params = load_params();
    
    if verbose {
        eprintln!("\nRunning with {} thread(s)!\n", &params.threads);
        // eprintln!("Params: {:?} ", &params);
        eprintln!("Processing bam:\n\t{}\n", &params.bam);
    }
    if params.method == "header" {
        let count_vec = countbamheader(&params);
        // let _none = writer_fn((&count_vec).to_vec(), "counts.txt.gz".to_string());
        let _none = writer_fn((&count_vec).to_vec(), params.outfile.to_string());
        eprintln!("\nDone!!\n\n");
        return;
    }
    // if params.method == "tag"{
    //     let count_vec = countbamheader(&params);
    //     let _none = writer_fn((&count_vec).to_vec(), "counts.txt.gz".to_string());
    //     eprintln!("\nDone!!\n\n");
    //     return;
    // }
}

// fn writer_fn (count_vec: Vec<Vec<u8>>, fname: String) -> Result<(), Box<dyn Error>> {
fn writer_fn (count_vec: Vec<String>, fname: String) -> Result<(), Box<dyn Error>> {
        let f = File::create(fname)?;
        let mut gz = GzBuilder::new()
                        .filename("counts_mm.txt.gz")
                        .write(f, Compression::default());
        for result in count_vec {
            gz.write_all(&result.as_bytes())?;
        }
        gz.finish()?;
        Ok(())
}




// fn remove_whitespace(s: &mut String) {
//     s.retain(|c| !c.is_whitespace());
// }

// fn countbamheader(params: &Params) -> Vec<Vec<u8>>{
fn countbamheader(params: &Params) -> Vec<String>{
    eprintln!("Processing using cb and umi in header");
    let split = params.stringsep.to_string();
    let ibam = &params.bam;
    let mut total: usize = 0;
    let mut err: usize = 0;
    let reader = bam::BamReader::from_path(ibam, params.threads).unwrap();
    let mut seqnames = Vec::new();
    let header = reader.header().clone();
    let hdata = header.reference_names();
    for seq in hdata {
        seqnames.push(seq);
        // eprintln!("{}", seq)
    }
    let mut data = Vec::new();
    
    for record in reader {
        let record = record.unwrap();
        let mut sample = "";
        total+=1;
        let readheader = match str::from_utf8(record.name()) {
            Ok(v) => v,
            Err(e) => panic!("\n\n*******Invalid UTF-8 sequence: {}*******\n\n", e),
        };
        let mut cb = match readheader.split(&split).nth(params.cbpos){
            Some(v) => v.to_string(),
            None => {
                err+=1;
                continue
            },
        };
        if cb.len() <= 16 {
            err+=1;
            continue
        }
        let mut umi = match readheader.split(&split).nth(params.umipos){
            Some(v) => v.to_string(),
            None => {
                err+=1;
                continue
            },
        };
        if params.sampleheader >= 0 {
            let sampleint = params.sampleheader as usize;
            // eprintln!("{}", sampleint);
            sample = match readheader.split(&split).nth(sampleint){
                Some(v) => v,
                None => {
                    err+=1;
                    continue
                },
            };
        }
        if umi.len() <= 10 {
            err+=1;
            continue
        }
        cb = cb.replace(&params.cbsep, "");
        umi = umi.replace(&params.umisep, "");
        // let ref_id = seqnames.iter().position(|&r| r == &record.ref_id());
        // let ref_id = seqnames[(record.ref_id() as usize)];
        let start = record.start();
        // let refid: usize = record.ref_id() as usize;
        let refid: String;
        if record.ref_id() == -1 {
            refid = "Unmapped".to_string();
        } else {
            let ind: usize = record.ref_id() as usize;
            refid = seqnames[ind].to_string();
        }
        let mut seq = "".to_string();
        if params.returnseq{
            if record.sequence().available() {
                // eprintln!("{:?}", record.sequence().to_vec());
                seq = String::from_utf8_lossy(&record.sequence().to_vec()).to_string();
            } else {
                err+=1;
                continue
            }
        }
        // eprintln!("{} {} {} {}", &cb, &umi, readheader, start);
        // data.push(format!("{} {} {} {}", &cb, &umi, record.ref_id(), start));
        // data.push(format!("{} {} {} {} {}", &cb, &umi, refid, start, record.mapq()));
        if params.sampleheader >= 0{
            if params.returnseq {
                data.push(format!("{} {} {} {} {} {} {}\n", sample, &cb, &umi, refid, start, record.mapq(), seq));
            } else {
                data.push(format!("{} {} {} {} {} {}\n", sample, &cb, &umi, refid, start, record.mapq()));
            }
        } else {
            if params.returnseq {
                data.push(format!("{} {} {} {} {} {}\n", &cb, &umi, refid, start, record.mapq(), seq));
            } else {
                data.push(format!("{} {} {} {} {}\n", &cb, &umi, refid, start, record.mapq()));
            }
            
        }

    }
    eprintln!("Found {} reads\n\tNumbers of errors: {}", total, err);
    data.sort();
    // let cdata = data;
    return data;
    // let mut out_vec = Vec::new();
    // let cdata = data.into_iter().dedup_with_count();
    // for (count, record) in cdata {
    //    let count_str = record + &" ".to_owned() + &(count.to_string()+&"\n".to_owned());
    //     out_vec.push(count_str.as_bytes().to_owned());
    // }
    // return out_vec;
}


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


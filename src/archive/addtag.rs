/**

time ~/develop/mutCaller/target/release/addtag --ibam Aligned.sortedByCoord.out.bam --obam Aligned.sortedByCoord.out.tagged.bam
time ~/develop/mutCaller/addTags.py -u 10 -c 16 Aligned.sortedByCoord.out.bam | samtools view -hbo kquant/Aligned.out.tagged.bam > Aligned.sortedByCoord.out.tagged.bam


**/

// extern crate clap;
// extern crate bam;
// extern crate simple_log;
// extern crate parasailors;

// use bam::RecordWriter;
// use clap::{App, load_yaml};
// use std::str;
// use simple_log::LogConfigBuilder;
// use std::io::{self};

#[macro_use]
extern crate simple_log;
extern crate clap;
extern crate bam;
// extern crate parasailors;

use simple_log::LogConfigBuilder;

use bam::RecordWriter;
use clap::{App, load_yaml};

use std::str;

// #[derive(Clone)]

struct Params {
    ibam: String,
    obam: String,
    split: String,
    joiner: String,
    threads: usize,
    // include_header: bool,
    cb_len: usize, 
    // umi_len: usize,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params_addtag.yml");
    let params = App::from_yaml(yaml).get_matches();
    let ibam = params.value_of("input_bam").unwrap();
    let obam = params.value_of("output_bam").unwrap_or("out.bam");
    let joiner = params.value_of("joiner").unwrap_or(":");
    // let include_header = params.value_of("no_header").is_some();
    // println!("include header: {:?}", params.value_of("no_header"));
    let split = params.value_of("split").unwrap_or("|BARCODE=");
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<usize>().unwrap() - 1;
    // let umi_len = params.value_of("umi_len").unwrap_or("10");
    // let umi_len = umi_len.to_string().parse::<usize>().unwrap() - 1;
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap();
    Params{
        ibam: ibam.to_string(),
        obam: obam.to_string(),
        threads: threads,
        split: split.to_string(),
        joiner: joiner.to_string(),
        // include_header: true,
        cb_len: cb_len,
        // umi_len: umi_len,
    }
}

fn main() {
    let params = load_params();
    let config = LogConfigBuilder::builder()
        .path("./addtags.log")
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("debug")
        .output_file()
        .build();

    let _ = simple_log::new(config);
    info!("starting!");
    eprintln!("Running with {} thread(s)!", params.threads);
    addtags(&params);
}


fn addtags(params: &Params) {
    // println!("include header: {:?}", *&params.include_header);
    // println!("threads: {:?}", *&params.threads);
    let (read_threads, write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as u16, ((*&params.threads/2) -1) as u16)
    } else {
        (0 as u16, 0 as u16)
    };
    let reader = bam::BamReader::from_path(params.ibam.to_string(), read_threads).unwrap();
    let mut writer = bam::BamWriter::build()
        .write_header(true)
        .additional_threads(write_threads)
        .from_path(&params.obam, reader.header().clone()).unwrap();

    for record in reader {
        let mut newrecord = record.unwrap();
        let seqname = match str::from_utf8(&newrecord.name()) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
        let cbumi= seqname.split(&params.split).nth(1).unwrap().to_string();
        let modified_name = seqname.replace(&params.split, &params.joiner);
        let (cb_umi_s1, cb_umi_s2) = cbumi.split_at(params.cb_len+1);
        // let mut newrecord = record.clone();
        newrecord.tags_mut().push_string(b"CB", &cb_umi_s1.as_bytes());
        newrecord.tags_mut().push_string(b"CR", &cb_umi_s1.as_bytes());
        newrecord.tags_mut().push_string(b"UB", &cb_umi_s2.as_bytes());
        newrecord.tags_mut().push_string(b"UR", &cb_umi_s2.as_bytes());
        newrecord.set_name(modified_name.bytes());
        writer.write(&newrecord).unwrap();
    }
}

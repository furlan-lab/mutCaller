/**

 **/


use vcf::{VCFReader, VCFHeaderFilterAlt, VCFError};
use flate2::read::MultiGzDecoder;
use flate2::GzBuilder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufReader, BufRead};
use crate::mutcaller::Variant;
use crate::utils::classify_variant;
use std::path::Path;
use std::error::Error;
use std::io::Write;
use rust_htslib::faidx::build;
use clap::{App, load_yaml};
use csv::ReaderBuilder;


pub fn faidx() {
    let yaml = load_yaml!("cli.yml");
    let matches = App::from_yaml(yaml).get_matches();
    let faidx_params = matches.subcommand_matches("FAIDX").unwrap();
    let fasta = faidx_params.value_of("genome").unwrap().to_string();
    eprintln!("Building fasta index for {}\n", fasta);
    let path = std::path::PathBuf::from(fasta);
    build(&path).expect("Failed to build fasta index");
}

pub fn check_variants_run() {
    let yaml = load_yaml!("cli.yml");
    let matches = App::from_yaml(yaml).get_matches();
    let cv_params = matches.subcommand_matches("CHECKVARIANTS").unwrap();
    let fasta = cv_params.value_of("genome").unwrap().to_string();
    let variantstring = cv_params.value_of("variants").unwrap().to_string();
    let mut classified_variants = Vec::new();
    let csvdata = read_csv(variantstring).unwrap();
    for variant in csvdata {
        let classified_variant = classify_variant(&variant);
        // eprintln!("Correctly parsed and classified variant: {}\n\n", classified_variant.as_ref().unwrap());
        classified_variants.push(classified_variant.unwrap());
    }
    let path = std::path::PathBuf::from(fasta);
    // build(&path).expect("Failed to build fasta index");
    let reader = rust_htslib::faidx::Reader::from_path(path).expect("Failed to open faidx");
    // assert_eq!(reader.seq_names(), Ok(vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()]));
    for variant in classified_variants {
        let seq = reader.fetch_seq_string(&variant.seq, variant.start.parse::<usize>().unwrap(), variant.start.parse::<usize>().unwrap()).unwrap();
        if seq == variant.ref_nt && variant.query_nt != variant.ref_nt {
            eprintln!("Correctly parsed, classified, and checked variant against reference genome: {}\n", variant);
        } else {
            eprintln!("Error in parsing, classifying, and checking variant against reference genome: {}\n", variant);
        }
        // eprintln!("Correctly parsed, classified, and checked variant against reference genome: {}\n", variant);
        // eprintln!("Found reference sequence: {} at position {} on {} \n", seq, variant.start, &variant.seq);
    }
}

fn read_csv(variantstring: String) -> Result<Vec<Variant>, Box<dyn Error>> {
    eprintln!("Opening variants file: {}\n", &variantstring);
    let file = File::open(variantstring).unwrap();
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


pub fn guess_vcf(file: &String) -> Result<bool, VCFError>{
    let path = Path::new(&file);
    let compression = {
        if path.extension().is_some() && path.extension().unwrap() == "gz" {
                true
            }else{
                false
        }
    };
    if compression {
        let mut reader = BufReader::new(MultiGzDecoder::new(File::open(&file)?));
        let mut first_line = String::new();
        let _ = reader.read_line(&mut first_line);
        first_line = String::from(first_line.trim());
        // eprintln!("{:?}", first_line);
        return Ok(first_line.contains("=VCF"));
    } else {
        let mut reader = BufReader::new(File::open(&file)?);
        let mut first_line = String::new();
        let _ = reader.read_line(&mut first_line);
        first_line = String::from(first_line.trim());
        // eprintln!("{:?}", first_line);
        return Ok(first_line.contains("=VCF"));
    }
}

pub fn guess_compression(file: &String) -> Result<bool, VCFError>{
    let path = Path::new(&file);
    if path.extension().is_some() && path.extension().unwrap() == "gz" {
            return Ok(true)
        }else{
            return Ok(false)
    };
}



pub fn read_vcf_compressed(file: &String, qual: &f64, verbose: &bool) -> Result<Vec<Variant>, VCFError> {
    let mut reader = VCFReader::new(BufReader::new(MultiGzDecoder::new(File::open(file)?)))?;

    // access FILTER contents
    assert_eq!(
        Some(VCFHeaderFilterAlt {
            id: b"PASS",
            description: b"All filters passed"
        }),
        reader.header().filter(b"PASS")
    );

    // prepare VCFRecord object
    let mut vcf_record = reader.empty_record();
    let mut count = 0;
    // read one record
    let mut data = Vec::new();
    while reader.next_record(& mut vcf_record).unwrap()  {
        count += 1;
        // todo count indels
        let alt = &vcf_record.alternative.clone().into_iter().flatten().collect::<Vec<u8>>();
        let rec = &vcf_record.reference;
        // let mut info_dat = Vec::new();
        // eprintln!("info: {:?}", info_dat);
        if vcf_record.qual.unwrap() > *qual {
            let mut vname = String::new();
            for (_l, v) in &vcf_record.info {
                let first_entry = v.into_iter().nth(0).unwrap();
                let split: Vec<_> = first_entry.split(|i| *i == 124).collect();

                vname = match split.len() {
                    4 => {
                        format!("{}_unknown", String::from_utf8_lossy(split[1]))
                    },
                    16 => {
                        format!("{}_{}_{}", String::from_utf8_lossy(split[3]), String::from_utf8_lossy(split[1]),String::from_utf8_lossy(split[9]))
                    },
                    _ => {
                        "name_not_determined".to_string()
                    }
                }
            }
            data.push(Variant{
                seq: String::from_utf8_lossy(&vcf_record.chromosome).to_string(),
                start: vcf_record.position.to_string(),
                ref_nt: String::from_utf8_lossy(rec).to_string(),
                query_nt: String::from_utf8_lossy(alt).to_string(),
                // ref_nt: String::from_utf8_lossy(rec).to_string().chars().nth(0).unwrap().to_string(),
                // query_nt: String::from_utf8_lossy(alt).to_string().chars().nth(0).unwrap().to_string(),
                name: vname.to_string(),
                class: None
            })
        } else {
            continue
        }
    }
    if *verbose {
        eprintln!("\tNumber of records in vcf file = {}\n\tAfter filtering, number of records = {}\n", count, data.len());
    }
    Ok(data)
    // data.push(Variant{seq: "chr12".to_string(),
    //             start: "12".to_string(),
    //             ref_nt: "A".chars().nth(0).unwrap(),
    //             query_nt: "C".chars().nth(0).unwrap(),
    //             name: "test".to_string()});
    // Ok(data)
}



pub fn read_vcf_uncompressed(file: &String, qual: &f64, verbose: &bool) -> Result<Vec<Variant>, VCFError> {
    let mut reader = VCFReader::new(BufReader::new(File::open(file)?))?;

    // access FILTER contents
    assert_eq!(
        Some(VCFHeaderFilterAlt {
            id: b"PASS",
            description: b"All filters passed"
        }),
        reader.header().filter(b"PASS")
    );

    // prepare VCFRecord object
    let mut vcf_record = reader.empty_record();
    let mut count = 0;
    // read one record
    let mut data = Vec::new();
    while reader.next_record(& mut vcf_record).unwrap()  {
        count += 1;
        // todo count indels
        let alt = &vcf_record.alternative.clone().into_iter().flatten().collect::<Vec<u8>>();
        let rec = &vcf_record.reference;
        // let mut info_dat = Vec::new();
        // eprintln!("info: {:?}", info_dat);
        if vcf_record.qual.unwrap() > *qual {
            let mut vname = String::new();
            for (_l, v) in &vcf_record.info {
                let first_entry = v.into_iter().nth(0).unwrap();
                let split: Vec<_> = first_entry.split(|i| *i == 124).collect();

                vname = match split.len() {
                    4 => {
                        format!("{}_unknown", String::from_utf8_lossy(split[1]))
                    },
                    16 => {
                        format!("{}_{}_{}", String::from_utf8_lossy(split[3]), String::from_utf8_lossy(split[1]),String::from_utf8_lossy(split[9]))
                    },
                    _ => {
                        "name_not_determined".to_string()
                    }
                }
            }
            data.push(Variant{
                seq: String::from_utf8_lossy(&vcf_record.chromosome).to_string(),
                start: vcf_record.position.to_string(),
                ref_nt: String::from_utf8_lossy(rec).to_string(),
                query_nt: String::from_utf8_lossy(alt).to_string(),
                // ref_nt: String::from_utf8_lossy(rec).to_string().chars().nth(0).unwrap().to_string(),
                // query_nt: String::from_utf8_lossy(alt).to_string().chars().nth(0).unwrap().to_string(),
                name: vname.to_string(),
                class: None
            })
        } else {
            continue
        }
    }
    if *verbose {
        eprintln!("\tNumber of records in vcf file = {}\n\tAfter filtering, number of records = {}\n", count, data.len());
    }
    Ok(data)
}


pub fn variants_writer_fn (file: String, variants: Vec<Variant>) -> Result<(), Box<dyn Error>> {
        let f = File::create(file.clone())?;
        let mut gz = GzBuilder::new()
                        .filename(file)
                        .write(f, Compression::default());
        gz.write_all(format!("seq\tstart\tref_nt\tquery_nt\tname\n").as_bytes())?;
        for variant in variants {
                // gz.write_all(&variant)?;
                gz.write(format!("{}\t{}\t{}\t{}\t{}\n", &variant.seq, &variant.start, &variant.ref_nt, &variant.query_nt, &variant.name).as_bytes())?;
        }
        gz.finish()?;
        Ok(())
}


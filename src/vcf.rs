/**

 **/


use vcf::{VCFReader, VCFHeaderFilterAlt, VCFError};
use flate2::read::MultiGzDecoder;
use flate2::GzBuilder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufReader, BufRead};
use crate::mutcaller::Variant;
use std::path::Path;
use std::error::Error;
use std::io::Write;


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


use csv::ReaderBuilder;
use faimm::IndexedFasta;
use flate2::read::{GzDecoder, MultiGzDecoder};
use flate2::{Compression, GzBuilder};
use serde::Deserialize;
use simplelog::info;
use std::error::Error;
use std::ffi::OsStr;
use std::path::PathBuf;
use std::fs::File;
use std::env;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use crate::mutcaller::{Paramsm, Variant, VariantClass};

pub fn get_current_working_dir() -> std::io::Result<PathBuf> {
    let cwd = env::current_dir();
    cwd
}


pub fn writer_fn(count_vec: Vec<Vec<Vec<u8>>>, params: &Paramsm) -> Result<(), Box<dyn Error>> {
    let counts_path = params.output_path.join("counts.txt.gz");
    let counts_file = counts_path.to_str().unwrap();
    info!("\n\n\tWriting counts to : '{}'\n", counts_file);
    if params.verbose {
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

pub fn check_variants(variant: &Variant, genome: &String, verbose: &bool) -> Variant {
    let start = variant.start.parse::<usize>().unwrap() - 1;
    let fa = IndexedFasta::from_file(genome).expect("Error opening fa");
    let chr_index = fa.fai().tid(&variant.seq).expect("Cannot find chr in index");
    let v = fa.view(chr_index, start, start + variant.ref_nt.len()).expect("Cannot get .fa view");
    assert_eq!(v.to_string(), variant.ref_nt, "Error in variants file; nt found at position: {} is {}; but nt provided was: {}", start + 1, v.to_string(), variant.ref_nt);
    info!("\tVariant ref_nt provided matches reference for the variant: '{}'", variant.name);
    if *verbose {
        eprintln!("Variant ref_nt matches reference for the variant: '{}'", variant.name);
    }
    variant.clone()
}

pub fn lines_from_file(filename: &str) -> Vec<String> {
    let path = Path::new(filename);
    let file = File::open(&path).expect(&format!("\n\n*******couldn't open {}*******\n\n", path.display()));
    let buf: Box<dyn BufRead> = if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    buf.lines().collect::<Result<_, _>>().expect("\n\n*******Could not parse line*******\n\n")
}


pub fn classify_variant(variant: &Variant) -> Result<Variant, Box<dyn Error>> {
    let mut classified_variant = variant.clone();
    classified_variant.class = if variant.ref_nt.len() == 1 && variant.query_nt.len() == 1 {
        Some(VariantClass::SNV)
    } else if variant.ref_nt.len() > 1 && variant.query_nt.len() > 1 {
        Some(VariantClass::MNV)
    } else if variant.ref_nt.len() == 1 && variant.query_nt.len() >= 1 {
        Some(VariantClass::Insertion)
    } else {
        Some(VariantClass::Deletion)
    };
    Ok(classified_variant)
}


pub fn read_csv(params: Option<&Paramsm>, file_i: Option<String>, verbose_i: Option<bool>) -> Result<Vec<Variant>, Box<dyn Error>> {
    let file = params.map_or_else(|| file_i.unwrap(), |p| p.variants.clone());
    let verbose = params.map_or_else(|| verbose_i.unwrap(), |p| p.verbose);

    #[derive(Deserialize)]
    struct PreVariant {
        seq: String,
        start: String,
        ref_nt: String,
        query_nt: String,
        name: String,
    }

    let path = Path::new(&file);
    let compression = path.extension().map_or(false, |ext| ext == "gz");

    let reader: Box<dyn BufRead> = if compression {
        Box::new(BufReader::new(MultiGzDecoder::new(File::open(&file)?)))
    } else {
        Box::new(BufReader::new(File::open(&file)?))
    };

    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(reader);

    let csvdata: Vec<Variant> = rdr.deserialize()
        .map(|result| {
            let prevariant: PreVariant = result?;
            Ok(Variant {
                seq: prevariant.seq,
                start: prevariant.start,
                ref_nt: prevariant.ref_nt,
                query_nt: prevariant.query_nt,
                name: prevariant.name,
                class: None,
            })
        })
        .collect::<Result<_, csv::Error>>()?;

    if verbose {
        eprintln!("Opening variants file: {}\n", file);
    }
    info!("Opening variants file: {}\n", file);
    Ok(csvdata)
}

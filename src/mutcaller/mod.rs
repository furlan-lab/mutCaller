extern crate simplelog;
extern crate csv;
extern crate clap;
extern crate bam;
extern crate serde;
extern crate fastq;
extern crate itertools;

// use bam::record::tags;
use clap::{App, load_yaml};
use itertools::Itertools;
use serde::Deserialize;
use simplelog::{Config, WriteLogger, CombinedLogger, LevelFilter, info, warn};
use std::error::Error;
use std::fmt;
use std::fs::{self, File};
use std::io::{Error as IoError, ErrorKind};
use std::path::{Path, PathBuf};
use std::str;
use std::time::Instant;
use crate::align::{align, Aligner, AlignerFlavor, test_progs};
use crate::countbam::Params;
use crate::utils::{read_csv, check_variants, classify_variant, writer_fn, get_current_working_dir};
use crate::vcf::{guess_vcf, guess_compression, read_vcf_compressed, read_vcf_uncompressed};
use crate::fastq::fastq;

#[derive(Deserialize, Debug, Clone)]
pub struct Variant {
    pub seq: String,
    pub start: String,
    pub ref_nt: String,
    pub query_nt: String,
    pub name: String,
    pub class: Option<VariantClass>,
}

#[derive(Deserialize, Debug, Clone, PartialEq)]
pub enum VariantClass {
    SNV,
    MNV,
    Insertion,
    Deletion,
    None,
}

impl fmt::Display for Variant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let class = self.class.as_ref().unwrap_or(&VariantClass::None);
        write!(f, "seq: {} start: {} ref_nt: {} query_nt: {} name: {} class: {:?}", self.seq, self.start, self.ref_nt, self.query_nt, self.name, class)
    }
}

#[derive(Debug, Clone)]
pub struct Paramsm {
    pub fastq1: String,
    pub fastq2: String,
    pub genome: String,
    pub bcs: String,
    pub umi_len: usize,
    pub cb_len: usize,
    pub threads: usize,
    pub aligner: Aligner,
    pub variants: String,
    pub read_len: usize,
    pub output_path: PathBuf,
    pub keep: bool,
    pub verbose: bool,
    pub vcf: bool,
    pub qual: f64,
}

fn load_params() -> Paramsm {
    let yaml = load_yaml!("../cli.yml");
    let matches = App::from_yaml(yaml).get_matches();
    let params = matches.subcommand_matches("UNALIGNED").unwrap();

    let fastq1 = params.value_of("fastq1").unwrap().to_string();
    let fastq2 = params.value_of("fastq2").unwrap().to_string();
    let output = params.value_of("output").unwrap_or("out");
    let genome = params.value_of("genome").unwrap().to_string();
    let bcs = params.value_of("barcodes_file").unwrap_or("/Users/sfurlan/develop/mutCaller/data/737K-august-2016.txt.gz").to_string();
    let threads = params.value_of("threads").unwrap_or("1").parse().unwrap();
    let umi_len = params.value_of("umi_len").unwrap_or("10").parse().unwrap();
    let cb_len = params.value_of("cb_len").unwrap_or("16").parse().unwrap();
    let read_len = params.value_of("read_len").unwrap_or("90").parse().unwrap();
    let aligner = params.value_of("aligner").unwrap_or("minimap2").to_string();
    let aligner_loc = params.value_of("aligner_loc").unwrap_or(&aligner).to_string();
    let add_aligner_args = params.value_of("add_aligner_args").map(|s| s.split_whitespace().map(str::to_string).collect());
    let qual = params.value_of("qual").unwrap_or("95.0").parse().unwrap();
    let verbose = params.is_present("verbose");
    let keep = params.is_present("keep_files");

    let aligner = Aligner::new(aligner, aligner_loc, add_aligner_args).unwrap();
    let variants = params.value_of("variants").unwrap().to_string();

    let params = Paramsm {
        fastq1,
        fastq2,
        genome,
        bcs,
        umi_len,
        cb_len,
        threads,
        aligner,
        variants,
        read_len,
        output_path: PathBuf::from(output),
        keep,
        verbose,
        vcf: false,
        qual,
    };
    let param_res = check_params(params);  
    let good_params = match param_res {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Error: {:?}", e);
            std::process::exit(1);
        }
    };
    
    return good_params
}


fn check_params(params: Paramsm) -> Result<Paramsm, Box<dyn Error>> {
    // Initialize logging
    let log_file_path = params.output_path.join("mutcaller.log");
    let log_file = log_file_path.to_str()
        .expect("Failed to convert log file path to string");
    
    CombinedLogger::init(vec![
        #[cfg(not(feature = "termcolor"))]
        WriteLogger::new(
            LevelFilter::Info,
            Config::default(),
            File::create(log_file).unwrap(),
        ),
    ]).unwrap();

    // CombinedLogger::init(vec![
    //     WriteLogger::new(LevelFilter::Info, Config::default(), File::create(log_file)?),
    // ])?;

    // Get current working directory
    let wdir = get_current_working_dir()?
        .to_str()
        .expect("Failed to convert working directory to string")
        .to_string();

    if params.verbose {
        eprintln!("\n\nCurrent working directory: '{:?}'", wdir);
    }

    // Resolve absolute output path
    let abs_outpath = if params.output_path.is_relative() {
        Path::new(&wdir).join(&params.output_path)
    } else {
        params.output_path.clone()
    };

    // Convert output path to string
    let abs_outpath_str = abs_outpath.to_str()
        .expect("Failed to convert absolute output path to string");

    // Check if output directory exists
    if abs_outpath.exists() {
        info!("\n\n\tFound existing output directory: '{}'\n", abs_outpath_str);
        warn!("\n\n\tExisting data in this folder could be lost!!!\n");
        if params.verbose {
            eprintln!("Found existing output directory: '{}'", abs_outpath_str);
            eprintln!("\tExisting data in this folder could be lost!!!");
        }
    } else {
        info!("\n\n\tCreating output directory: '{}'\n", abs_outpath_str);
        if params.verbose {
            eprintln!("Creating output directory: '{}'", abs_outpath_str);
        }
        fs::create_dir(&abs_outpath)?;
    }

    // Guess VCF format
    let is_vcf = guess_vcf(&params.variants)?;

    // Check that genome file exists
    let genome_path = Path::new(&params.genome);
    if !genome_path.exists() {
        let genome_path_str = genome_path.to_str()
            .expect("Failed to convert genome path to string");
        let msg = format!("Genome file '{}' does not exist", genome_path_str);
        eprint!("Error: {}\n", msg);
        return Err(IoError::new(ErrorKind::NotFound, msg).into());
    }

    Ok(Paramsm {
        vcf: is_vcf,
        ..params
    })
}


pub fn mutcaller_run() {
    // eprintln!("Running mutcaller");
    let start = Instant::now();
    let params = load_params();
    info!("\n\n\tParsing Parameters!\n");
    if params.verbose {
        eprintln!("\n\nParsing Parameters!\n");
    }
    info!("\n\n\tChecking programs and parsing variants!\n");
    if params.verbose {
        eprintln!("\n\nChecking programs and parsing variants!\n");
    }

    let _prog_test_res = test_progs();
    let csvdata = if params.vcf {
        let is_compressed = guess_compression(&params.variants).unwrap();
        if is_compressed {
            read_vcf_compressed(&params.variants, &params.qual, &params.verbose)
        } else {
            read_vcf_uncompressed(&params.variants, &params.qual, &params.verbose)
        }
    } else {
        Ok(read_csv(Some(&params), None, None).unwrap())
    };

    let checked_and_classified_variants: Vec<_> = csvdata.as_ref().unwrap().iter()
        .map(|variant| {
            let checked_variant = if params.aligner.flavor == AlignerFlavor::Minimap2 {
                check_variants(variant, &params.genome, &params.verbose)
            } else {
                variant.clone()
            };
            let classified_variant = classify_variant(&checked_variant).unwrap();
            if params.verbose {
                eprintln!("Correctly parsed and classified variant: {}\n\n", classified_variant);
            }
            info!("\tCorrectly parsed and classified variant: {}\n\n", classified_variant);
            classified_variant
        })
        .collect();

    info!("\n\n\tRunning with {} thread(s)!\n", params.threads);
    if params.verbose {
        eprintln!("\n\nRunning with {} thread(s)!\n", params.threads);
    }
    let _fqr = fastq(&params);
    info!("done!");
    let _ar = align(&params);

    let count_vec: Vec<_> = checked_and_classified_variants.into_iter()
        .filter_map(|variant| Some(count_variants_helper(Some(&params), None, variant.clone()).unwrap_or_else(|| fallback_count_variants(&variant))))
        .collect();

    let _none = writer_fn(count_vec, &params);
    let duration = start.elapsed();
    info!("\n\n\tDone!!\n");
    info!("\n\n\tTime elapsed is: {:?}\n", duration);
    if params.verbose {
        eprintln!("\n\nDone!!");
        eprintln!("\n\nTime elapsed is: {:?}", duration);
    }
}

fn process_variant(ref_id: u32, start: u32) -> bam::Region {
    bam::Region::new(ref_id, start - 1, start - 1)
}

pub fn count_variants_helper(paramsm: Option<&Paramsm>, params: Option<&Params>, variant: Variant) -> Option<Vec<Vec<u8>>> {
    let result = match variant.class.clone()? {
        VariantClass::SNV => {
            let verbose = paramsm.map_or_else(|| params.unwrap().verbose, |p| p.verbose);
            let cb_len = paramsm.map_or_else(|| None, |p| Some(p.cb_len));
            let ibam_temp = paramsm.map_or_else(|| Some(params.unwrap().bam.to_string()), |p| Some(p.output_path.join("Aligned.sortedByCoord.out.bam").to_str().unwrap().to_string()));
            let tags = params.map_or_else(|| None, |p| Some((p.cb_tag.as_bytes(), p.umi_tag.as_bytes())));
            if tags.is_none() {
                count_variants_snv(ibam_temp.as_deref(), verbose, cb_len, variant.clone(), Some("|BARCODE=".to_string()), None)
            } else {
                count_variants_snv(ibam_temp.as_deref(), verbose, cb_len, variant.clone(), None, tags)
            }

        }
        VariantClass::Deletion | VariantClass::Insertion => {
            let verbose = paramsm.map_or_else(|| params.unwrap().verbose, |p| p.verbose);
            let cb_len = paramsm.map_or_else(|| None, |p| Some(p.cb_len));
            let ibam_temp = paramsm.map_or_else(|| Some(params.unwrap().bam.to_string()), |p| Some(p.output_path.join("Aligned.sortedByCoord.out.bam").to_str().unwrap().to_string()));
            let tags = params.map_or_else(|| None, |p| Some((p.cb_tag.as_bytes(), p.umi_tag.as_bytes())));
            if tags.is_none() {
                count_variants_indel(ibam_temp.as_deref(), verbose, cb_len, variant.clone(), Some("|BARCODE=".to_string()), None)
            } else {
                count_variants_indel(ibam_temp.as_deref(), verbose, cb_len, variant.clone(), None, tags)
            }
        }
        _ => {
            warn!("\n\n\tVariant type {:?} not currently supported", variant.class.as_ref().unwrap());
            if paramsm.unwrap().verbose {
                eprintln!("\n\n\tVariant type {:?} not currently supported", variant.class.clone().unwrap());
            }
            None
        }
    };

    result.or_else(|| Some(fallback_count_variants(&variant)))
}

fn fallback_count_variants(_variant: &Variant) -> Vec<Vec<u8>> {
    vec![b"fallback data".to_vec()]
}

fn string_pop(slice: &[u8]) -> &[u8; 2] {
    slice.try_into().expect("slice with incorrect length")
}

fn get_cb(split: Option<String>, cb_len: Option<usize>, readname: String, record: &bam::Record, tags: Option<(&[u8], &[u8])>) -> Result<(String, String), IoError> {
    // let barcode_error = IoError::new(ErrorKind::Other, "Cell barcode / UMI not found");
    // eprintln!("split: {:?}, cb_len: {:?}, readname: {:?}, record: {:?}, tags: {:?}", split, cb_len, readname, record, tags);
    if let Some(split_str) = split {
        let cbumi = readname.split(&split_str).nth(1).ok_or_else(|| IoError::new(ErrorKind::Other, "Cell barcode / UMI not found"))?;
        let (cb, umi) = cbumi.split_at(cb_len.unwrap() + 1);
        Ok((cb.to_string(), umi.to_string()))
    } else {
        let cb_tag_b = string_pop(tags.unwrap().0);
        let umi_tag_b = string_pop(tags.unwrap().1);
        let cb = record.tags().get(cb_tag_b).and_then(|tag| match tag {
            bam::record::tags::TagValue::String(cba, _) => Some(str::from_utf8(&cba).unwrap().to_string()),
            _ => None,
        }).ok_or_else(|| IoError::new(ErrorKind::Other, "Cell barcode / UMI not found"))?;
        let umi = record.tags().get(umi_tag_b).and_then(|tag| match tag {
            bam::record::tags::TagValue::String(uma, _) => Some(str::from_utf8(&uma).unwrap().to_string()),
            _ => None,
        }).ok_or_else(|| IoError::new(ErrorKind::Other, "Cell barcode / UMI not found"))?;
        Ok((cb, umi))
    }
}

fn get_header_seqs(header: bam::Header) -> Vec<String> {
    header.reference_names().iter().map(|seq| seq.to_string()).collect()
}

#[derive(Debug, PartialEq, Clone)]
enum MatchType {
    Ref,
    Query,
    Other,
}

#[derive(PartialEq, Clone)]
struct SequenceMatch {
    still_to_check: usize,
    final_result: Option<MatchType>,
}

fn count_variants_snv(
    ibam: Option<&str>, 
    verbose: bool, 
    cb_len: Option<usize>, 
    variant: Variant, 
    split: Option<String>, 
    tags: Option<(&[u8], &[u8])>
) -> Option<Vec<Vec<u8>>> {
    let debug = false;
    let ibam = ibam?;
    let seqname = variant.seq;
    let start = variant.start.parse::<u32>().ok()?;
    let vname = variant.name;
    let query_nt = variant.query_nt.chars().next()?;
    let us_ref_nt = variant.ref_nt.chars().next()?;

    let mut reader = bam::IndexedReader::build().from_path(&ibam).ok()?;
    let seqnames = get_header_seqs(reader.header().clone());
    let ref_id = seqnames.iter().position(|r| r == &seqname)?;
    let region = process_variant(ref_id as u32, start);

    let mut data = Vec::new();
    let mut total = 0;
    let mut err = 0;
    let mut query = 0;
    let mut reference = 0;
    let mut other = 0;

    for record in reader.fetch(&region).ok()? {
        total += 1;
        let readname = str::from_utf8(record.as_ref().unwrap().name()).ok()?;
        if let Ok((cb, umi)) = get_cb(split.clone(), cb_len, readname.to_string(), record.as_ref().unwrap(), tags) {
            for entry in record.as_ref().unwrap().alignment_entries().ok()? {
                // eprintln!("entry: {:?}", entry);
                if let Some((ref_pos, _ref_nt)) = entry.ref_pos_nt() {
                    if ref_pos == start - 1 {
                        if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
                            let match_type = if record_nt as char == us_ref_nt && !entry.is_insertion() && !entry.is_deletion() {
                                reference += 1;
                                if debug{
                                    eprintln!("reference nt = {}\n record_nt = {}\n ergo: reference", _ref_nt as char, record_nt as char);
                                };
                                MatchType::Ref
                            } else if record_nt as char == query_nt && !entry.is_insertion() && !entry.is_deletion() {
                                query += 1;
                                if debug{
                                    eprintln!("reference nt = {}\n record_nt = {}\n ergo: query", _ref_nt as char, record_nt as char);
                                };
                                MatchType::Query
                            } else if record_nt == b'N' {
                                err += 1;
                                continue;
                            } else {
                                other += 1;
                                if debug{
                                    eprintln!("reference nt = {}\n record_nt = {}\n ergo: Other", _ref_nt as char, record_nt as char);
                                };
                                MatchType::Other
                            };
                            
                            data.push(format!("{} {} {} {} {} {:?}", cb, umi, seqname, start, vname, match_type));
                        }
                    }
                }
            }
        } else {
            eprintln!("Error: Cell barcode / UMI not found");
            err += 1;
        }
    }

    log_counts(total, vname, err, query, reference, other, verbose);

    data.sort();
    let mut out_vec = Vec::new();
    for (count, record) in data.into_iter().dedup_with_count() {
        out_vec.push(format!("{} {}\n", record, count).as_bytes().to_owned());
    }

    Some(out_vec)
}


fn count_variants_indel(
    ibam: Option<&str>, 
    verbose: bool, 
    cb_len: Option<usize>, 
    variant: Variant, 
    split: Option<String>, 
    tags: Option<(&[u8], &[u8])>
) -> Option<Vec<Vec<u8>>> {
    let ibam = ibam?;
    let seqname = variant.seq.clone();
    let start = variant.start.parse::<u32>().ok()?;
    let vname = variant.name.clone();
    let class = variant.class.clone()?;

    let mut reader = bam::IndexedReader::build().from_path(&ibam).ok()?;
    let seqnames = get_header_seqs(reader.header().clone());
    let ref_id = seqnames.iter().position(|r| r == &seqname)?;
    let region = process_variant(ref_id as u32, start);

    let mut data = Vec::new();
    let mut total = 0;
    let mut err = 0;
    let mut query = 0;
    let mut reference = 0;
    let mut other = 0;

    for record in reader.fetch(&region).ok()? {
        total += 1;
        let readname = str::from_utf8(record.as_ref().unwrap().name()).ok()?;
        if let Ok((cb, umi)) = get_cb(split.clone(), cb_len, readname.to_string(), record.as_ref().unwrap(), tags) {
            let mut ins_result = SequenceMatch { still_to_check: variant.query_nt.len(), final_result: None };
            let mut del_result = SequenceMatch { still_to_check: variant.ref_nt.len(), final_result: None };
            for entry in record.as_ref().unwrap().alignment_entries().ok()? {
                match class {
                    VariantClass::Insertion => handle_insertion(entry, start, &variant, &mut ins_result, &mut query, &mut reference, &mut other),
                    VariantClass::Deletion => handle_deletion(entry, start, &variant, &mut del_result, &mut query, &mut reference, &mut other),
                    _ => continue,
                }

                if let Some(result) = ins_result.final_result.clone().or(del_result.final_result.clone()) {
                    data.push(format!("{} {} {} {} {} {:?}", cb, umi, seqname, start, vname, result));
                    break;
                }
            }
        } else {
            err += 1;
        }
    }

    log_counts(total, vname, err, query, reference, other, verbose);

    data.sort();
    let mut out_vec = Vec::new();
    for (count, record) in data.into_iter().dedup_with_count() {
        out_vec.push(format!("{} {}\n", record, count).as_bytes().to_owned());
    }

    Some(out_vec)
}

fn handle_insertion(
    entry: bam::record::AlignmentEntry,
    start: u32,
    variant: &Variant,
    result: &mut SequenceMatch,
    query: &mut usize,
    reference: &mut usize,
    other: &mut usize
) {
    if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
        if ref_pos >= start - 1 {
            if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
                update_result(record_nt, ref_nt, &variant.query_nt, result, query, reference, other, entry.is_insertion());
            }
        }
    }
}

fn handle_deletion(
    entry: bam::record::AlignmentEntry,
    start: u32,
    variant: &Variant,
    result: &mut SequenceMatch,
    query: &mut usize,
    reference: &mut usize,
    other: &mut usize
) {
    if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
        if ref_pos >= start - 1 {
            if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
                update_result(record_nt, ref_nt, &variant.ref_nt, result, query, reference, other, entry.is_deletion());
            }
        }
    }
}

fn update_result(
    record_nt: u8,
    ref_nt: u8,
    variant_nt: &str,
    result: &mut SequenceMatch,
    query: &mut usize,
    reference: &mut usize,
    other: &mut usize,
    is_indel: bool
) {
    if variant_nt.len() == result.still_to_check {
        result.still_to_check -= 1;
        if record_nt as char != ref_nt as char {
            result.final_result = Some(MatchType::Other);
            *other += 1;
        }
    } else if record_nt as char == variant_nt.chars().nth(variant_nt.len() - result.still_to_check).unwrap() && is_indel {
        result.still_to_check -= 1;
        if result.still_to_check == 0 {
            result.final_result = Some(MatchType::Query);
            *query += 1;
        }
    } else if record_nt as char == ref_nt as char && !is_indel {
        result.still_to_check -= 1;
        if result.still_to_check == 0 {
            result.final_result = Some(MatchType::Ref);
            *reference += 1;
        }
    } else {
        result.still_to_check -= 1;
        if result.still_to_check == 0 {
            result.final_result = Some(MatchType::Other);
            *other += 1;
        }
    }
}

fn log_counts(total: usize, vname: String, err: usize, query: usize, reference: usize, other: usize, verbose: bool) {
    info!("\n\n\tFound {} reads spanning variant: {}!\n\tNumbers of errors: {}", total, vname, err);
    if verbose {
        eprintln!("Found {} reads spanning variant: {}!\n\tNumbers of errors: {}", total, vname, err);
    }
    info!("\t\t\tQuery: {}\n\t\t\tReference: {}\n\t\t\tOther: {}", query, reference, other);
    if verbose {
        eprintln!("\t\tQuery: {}\n\t\tReference: {}\n\t\tOther: {}", query, reference, other);
    }
}

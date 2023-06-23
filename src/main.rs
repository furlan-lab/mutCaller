
use clap::{App, load_yaml};


pub mod mutcaller;
pub mod countbam;
pub mod vcf;

use crate::mutcaller::mutcaller_run;
use crate::countbam::{countbam_run};
use crate::vcf::{read_vcf_compressed, read_vcf_uncompressed, guess_vcf, guess_compression, variants_writer_fn};
use crate::mutcaller::read_csv;
use clap::ArgMatches;
// use itertools::Itertools;

fn main() {
    let yaml = load_yaml!("cli.yml");
    let params = App::from_yaml(yaml).get_matches();
    if let Some(params) = params.subcommand_matches("ARGPARSE") {
        argparse_run(&params);
    }
    if let Some(_params) = params.subcommand_matches("UNALIGNED") {
    	mutcaller_run();
	}
	if let Some(_params) = params.subcommand_matches("ALIGNED") {
    	countbam_run()
	}
    if let Some(params) = params.subcommand_matches("VARIANTS") {
        let qual = params.value_of("qual").unwrap_or("95.0");
        let mut verbose = true;
        if params.is_present("quiet") {
                verbose = false
        };
        let mut make_variants = "false";
        if params.is_present("make_variants") {
                make_variants = params.value_of("make_variants").unwrap();
        };
        let qual_p = qual.parse::<f64>();
        let vcf_file = params.value_of("variants").unwrap_or("/Users/sfurlan/develop/mutCaller/tests/var.vcf.gz").to_string();
        if verbose {
            eprintln!("Reading variants file: {}\n", &vcf_file);
        }
        let is_vcf = guess_vcf(&vcf_file);
        let is_compressed = guess_compression(&vcf_file);
        let data = {
            if is_vcf.unwrap(){
                if is_compressed.unwrap(){
                    read_vcf_compressed(&vcf_file.to_string(), &qual_p.unwrap(), &verbose).unwrap() 
               } else {
                    read_vcf_uncompressed(&vcf_file.to_string(), &qual_p.unwrap(), &verbose).unwrap() 
               }
            } else {
                read_csv(None, Some(vcf_file), Some(true)).unwrap()
            }
        };
        eprintln!("Records:\n");
        for variant in &data {
            eprintln!("{}", variant);
        }
        if make_variants != "false" {
            let res = variants_writer_fn(make_variants.to_string(), data);
            if verbose && res.unwrap()== () {
                eprintln!("Created variants file: {}\n", make_variants);
            }
        }   
    }

}



fn argparse_run(params: &ArgMatches){
    let arg1 = if params.is_present("arg1"){
        params.value_of("arg1").unwrap()
    } else {
        "arg1 not given"
    };
    let arg2 = if params.is_present("arg2"){
        params.value_of("arg2").unwrap()
    } else {
        "arg2 not given"
    };

    eprintln!("Argument 1: {:?}; Argument 2: {:?}", arg1, arg2);
    eprintln!("Argument 1 split by space");
    for byte in arg1.split_whitespace(){
        println!("{}", byte);
    }
    eprintln!("Argument 2 split by space");
    for byte in arg2.split_whitespace(){
        println!("{}", byte);
    }
}
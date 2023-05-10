
// #[macro_use] extern crate arrayref;

// use simple_log::LogConfigBuilder;
// use simple_log::info;
use clap::{App, load_yaml};

pub mod mutcaller;
pub mod countbam;
use crate::mutcaller::mutcaller_run;
use crate::countbam::countbam_run;

fn main() {

    let yaml = load_yaml!("cli.yml");
    let params = App::from_yaml(yaml).get_matches();
    if let Some(params) = params.subcommand_matches("UNALIGNED") {
    	eprintln!("Running unaligned with params: {:?}", params);
    	mutcaller_run();
	}
	if let Some(params) = params.subcommand_matches("ALIGNED") {
    	eprintln!("Running aligned with params: {:?}", params);
    	countbam_run()
	}

}

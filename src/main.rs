// #[macro_use] extern crate log;
// use simple_log::LogConfigBuilder;
// use simple_log::info;
use clap::{App, load_yaml};
// use simple_log::LogConfigBuilder;

pub mod mutcaller;
pub mod countbam;
use crate::mutcaller::mutcaller_run;
use crate::countbam::{countbam_run};

fn main() {
	// let wdpb= get_current_working_dir().unwrap();
    // let wdir = wdpb.to_str().unwrap();
    let yaml = load_yaml!("cli.yml");
    let params = App::from_yaml(yaml).get_matches();
    // let checked_params = countbam::load_params();
    // let config = LogConfigBuilder::builder()
    //     .path(checked_params.output_path.join("mutcaller.log").to_str().unwrap())
    //     .size(1 * 100)
    //     .roll_count(10)
    //     .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
    //     .level("info")
    //     .output_file()
    //     .build();
    // let _ = simple_log::new(config);
    // if checked_params.verbose {
    //     eprintln!("Writing log file: '{}'", &checked_params.output_path.join("mutcaller.log").to_str().unwrap());
    // }

    // if checked_params.verbose {
    //     eprintln!("\n\nCurrent working directory: '{}'", wdir);
    // }
    if let Some(_params) = params.subcommand_matches("UNALIGNED") {
    	// eprintln!("Running unaligned with params: {:?}", params);
    	mutcaller_run();
	}
	if let Some(_params) = params.subcommand_matches("ALIGNED") {
    	// eprintln!("Running aligned with params: {:?}", params);
    	countbam_run()
	}

}

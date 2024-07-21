extern crate simplelog;
extern crate csv;
extern crate clap;
extern crate bam;
extern crate serde;
extern crate fastq;
extern crate itertools;
extern crate rayon;



use std::fs;
use std::io;
use std::error::Error;
use std::process::Command;
use std::fs::File;
use std::process::Stdio;
use serde::Deserialize;
use crate::mutcaller::Paramsm;
use std::io::{Error as IoError, ErrorKind};
use std::path::Path;

#[cfg(not(feature = "paris"))]
use log::*;

#[derive(Deserialize, Clone, Debug, PartialEq)]
pub struct Aligner {
    pub flavor: AlignerFlavor,
    pub loc: String,
    pub args: Option<Vec<String>>,
}

#[derive(Deserialize, Clone, Debug, PartialEq)]
pub enum AlignerFlavor {
    Minimap2,
    STAR,
    BWA,
}



impl Aligner {
    pub fn new(aligner: String, aligner_loc: String, aligner_args: Option<Vec<String>>) -> Result<Aligner, io::Error> {
        let aligner_error = IoError::new(ErrorKind::Other, "Aligner not configured");
        Command::new(&aligner_loc)
            .arg("--help")
            .stderr(Stdio::piped())
            .stdout(Stdio::piped())
            .output()
            .expect("\n\n*******Failed to execute aligner*******\n\n");
            // .map_err(|_| aligner_error.clone())?;

        match aligner.as_str() {
            "minimap2" => Ok(Aligner { flavor: AlignerFlavor::Minimap2, loc: aligner_loc, args: aligner_args }),
            "bwa" => Ok(Aligner { flavor: AlignerFlavor::BWA, loc: aligner_loc, args: aligner_args }),
            "STAR" => Ok(Aligner { flavor: AlignerFlavor::STAR, loc: aligner_loc, args: aligner_args }),
            _ => Err(aligner_error),
        }
    }
}

pub fn test_progs() -> Result<(), Box<dyn Error>> {
    Command::new("samtools")
        .arg("-h")
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .output()
        .expect("\n\n*******Failed to execute samtools*******\n\n");
    Ok(())
}

// pub fn align(params: &Paramsm) -> Result<(), Box<dyn Error>> {
//     let align_sam = params.output_path.join("Aligned.sam");
//     let align_sorted_sam = params.output_path.join("Aligned.sorted.sam");
//     let align_sorted_bam = params.output_path.join("Aligned.sortedByCoord.out.bam");
//     let fastq = params.output_path.join("mutcaller_R1.fq.gz");
//     let outfolder = params.output_path.join("");
//     let zcat_cmd = if cfg!(target_os = "macos") { "zcat <" } else { "zcat" };

//     match params.aligner.flavor {
//         AlignerFlavor::Minimap2 => {
//             let args = params.aligner.args.clone();
//             if params.verbose {
//                 eprintln!("Aligning reads using minimap2");
//             }
//             info!("Aligning reads using minimap2");
//             let output = Command::new(&params.aligner.loc)
//                 .args(["--MD", "-Y", "-a", &params.genome, "-t", &params.threads.to_string(), "-o", align_sam.to_str().unwrap()])
//                 .args(args.unwrap_or_default())
//                 .arg(fastq.to_str().unwrap())
//                 .stderr(Stdio::piped())
//                 .stdout(Stdio::piped())
//                 .output()
//                 .expect("\n\n*******Failed to execute minimap2*******\n\n");
//             if params.verbose {
//                 eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//                 eprintln!("Minimap2 complete; Running samtools sort");
//             }
//             info!("{}", String::from_utf8_lossy(&output.stderr));
//             info!("Minimap2 complete; Running samtools sort");
//         }
//         AlignerFlavor::STAR => {
//             let args = params.aligner.args.clone();
//             if params.verbose {
//                 eprintln!("Aligning reads using STAR");
//             }
//             info!("Aligning reads using STAR");
//             let output = Command::new(&params.aligner.loc)
//                 .arg("--outFileNamePrefix").arg(outfolder.to_str().unwrap())
//                 .arg("--genomeDir").arg(&params.genome)
//                 .arg("--readFilesIn").arg(fastq.to_str().unwrap())
//                 .arg("--readNameSeparator").arg("space")
//                 .arg("--runThreadN").arg(params.threads.to_string())
//                 .arg("--outSAMunmapped").arg("Within").arg("KeepPairs")
//                 .arg("--outSAMtype").arg("BAM").arg("SortedByCoordinate")
//                 .arg("--outSAMattributes").arg("All")
//                 .arg("--readFilesCommand").arg(zcat_cmd)
//                 .args(args.unwrap_or_default())
//                 .stderr(Stdio::piped())
//                 .stdout(Stdio::piped())
//                 .output()
//                 .expect("\n\n*******Failed to execute STAR*******\n\n");
//             if params.verbose {
//                 eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//                 eprintln!("STAR complete; running samtools index");
//             }
//         }
//         AlignerFlavor::BWA => {
//             let stdout = File::create(align_sam.to_str().unwrap()).expect("failed to open output file for bwa mem");
//             let args = params.aligner.args.clone();
//             if params.verbose {
//                 eprintln!("Aligning reads using bwa mem");
//             }
//             info!("Aligning reads using bwa mem");
//             let output = Command::new(&params.aligner.loc)
//                 .arg("mem")
//                 .arg("-t").arg(params.threads.to_string())
//                 .args(args.unwrap_or_default())
//                 .arg(&params.genome)
//                 .arg(fastq.to_str().unwrap())
//                 .stderr(Stdio::piped())
//                 .stdout(stdout)
//                 .output()
//                 .expect("\n\n*******Failed to execute bwa*******\n\n");
//             if params.verbose {
//                 eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//                 eprintln!("bwa mem complete; Running samtools sort");
//             }
//             info!("{}", String::from_utf8_lossy(&output.stderr));
//             info!("bwa mem complete; Running samtools sort");
//         }
//     }

//     let output = Command::new("samtools")
//         .arg("sort")
//         .arg("-@").arg(params.threads.to_string())
//         .arg("-o").arg(align_sorted_sam.to_str().unwrap())
//         .arg(align_sam.to_str().unwrap())
//         .stderr(Stdio::piped())
//         .stdout(Stdio::piped())
//         .output()
//         .expect("\n\n*******Failed to execute samtools view*******\n\n");
//     if params.verbose {
//         eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//         eprintln!("Samtools sort complete; Running samtools view");
//     }
//     info!("{}", String::from_utf8_lossy(&output.stderr));
//     info!("Samtools sort complete; Running samtools view");

//     let output = Command::new("samtools")
//         .arg("view")
//         .arg("-b")
//         .arg("-@").arg(params.threads.to_string())
//         .arg("-o").arg(align_sorted_bam.to_str().unwrap())
//         .arg(align_sorted_sam.to_str().unwrap())
//         .stderr(Stdio::piped())
//         .stdout(Stdio::piped())
//         .output()
//         .expect("\n\n*******Failed to execute samtools sort*******\n\n");
//     if params.verbose {
//         eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//         eprintln!("Samtools view complete; Running samtools index");
//     }
//     info!("{}", String::from_utf8_lossy(&output.stderr));
//     info!("Samtools view complete; Running samtools index");

//     let output = Command::new("samtools")
//         .arg("index")
//         .arg("-@").arg(params.threads.to_string())
//         .arg(align_sorted_bam.to_str().unwrap())
//         .stdout(Stdio::piped())
//         .stderr(Stdio::piped())
//         .output()
//         .expect("\n\n*******Failed to execute samtools index*******\n\n");
//     if params.verbose {
//         eprintln!("{}", String::from_utf8_lossy(&output.stderr));
//     }
//     info!("{}", String::from_utf8_lossy(&output.stderr));

//     if !params.keep {
//         match params.aligner.flavor {
//             AlignerFlavor::Minimap2 => {
//                 fs::remove_file(align_sorted_sam.to_str().unwrap())?;
//                 fs::remove_file(align_sam.to_str().unwrap())?;
//                 fs::remove_file(fastq.to_str().unwrap())?;
//             }
//             AlignerFlavor::STAR => {
//                 fs::remove_file(fastq.to_str().unwrap())?;
//                 let possible_tmp = params.output_path.join("_STARtmp");
//                 if possible_tmp.exists() {
//                     fs::remove_dir_all(possible_tmp)?;
//                 }
//             }
//             _ => (),
//         }
//     }

//     Ok(())
// }


pub fn align(params: &Paramsm) -> Result<(), Box<dyn Error>> {
    let align_sam = params.output_path.join("Aligned.sam");
    let align_sorted_sam = params.output_path.join("Aligned.sorted.sam");
    let align_sorted_bam = params.output_path.join("Aligned.sortedByCoord.out.bam");
    let fastq = params.output_path.join("mutcaller_R1.fq.gz");
    let outfolder = params.output_path.join("");

    let zcat_cmd = if cfg!(target_os = "macos") { "zcat <" } else { "zcat" };

    match params.aligner.flavor {
        AlignerFlavor::Minimap2 => {
            run_minimap2(&params, &align_sam, &align_sorted_sam, &align_sorted_bam, &fastq)?;
        }
        AlignerFlavor::STAR => {
            run_star(&params, &outfolder, &fastq, zcat_cmd)?;
        }
        AlignerFlavor::BWA => {
            run_bwa(&params, &align_sam, &align_sorted_sam, &align_sorted_bam, &fastq)?;
        }
    }

    index_bam(&params, &align_sorted_bam)?;

    if !params.keep {
        cleanup(&params, &align_sam, &align_sorted_sam, &fastq)?;
    }

    Ok(())
}

fn run_minimap2(params: &Paramsm, align_sam: &Path, align_sorted_sam: &Path, align_sorted_bam: &Path, fastq: &Path) -> Result<(), Box<dyn Error>> {
    let command = &params.aligner.loc;
    let args = params.aligner.args.clone().unwrap_or_default();
    log_info(params, "Aligning reads using minimap2");

    let output = Command::new(command)
        .args(&["--MD", "-Y", "-a", &params.genome.to_string(), "-t", &params.threads.to_string(), "-o", align_sam.to_str().unwrap()])
        .args(&args)
        .arg(fastq.to_str().unwrap())
        .output()
        .expect("Failed to execute minimap2");

    log_output(params, &output.stderr, "Minimap2 complete; Running samtools sort");

    run_samtools_sort(params, align_sam, align_sorted_sam)?;
    run_samtools_view(params, align_sorted_sam, align_sorted_bam)
}

fn run_star(params: &Paramsm, outfolder: &Path, fastq: &Path, zcat_cmd: &str) -> Result<(), Box<dyn Error>> {
    let command = &params.aligner.loc;
    let args = params.aligner.args.clone().unwrap_or_default();
    log_info(params, "Aligning reads using STAR");

    let output = Command::new(command)
        .args(&["--outFileNamePrefix", outfolder.to_str().unwrap(), "--genomeDir", &params.genome.to_string(), "--readFilesIn", fastq.to_str().unwrap(), "--readNameSeparator", "space", "--runThreadN", &params.threads.to_string(), "--outSAMunmapped", "Within", "KeepPairs", "--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMattributes", "All", "--readFilesCommand", zcat_cmd])
        .args(&args)
        .output()
        .expect("Failed to execute STAR");

    log_output(params, &output.stderr, "STAR complete; running samtools index");
    Ok(())
}

fn run_bwa(params: &Paramsm, align_sam: &Path, align_sorted_sam: &Path, align_sorted_bam: &Path, fastq: &Path) -> Result<(), Box<dyn Error>> {
    let command = &params.aligner.loc;
    let args = params.aligner.args.clone().unwrap_or_default();
    log_info(params, "Aligning reads using bwa mem");

    let stdout = File::create(align_sam).expect("Failed to open output file for bwa mem");

    let output = Command::new(command)
        .args(&["mem", "-t", &params.threads.to_string(), &params.genome.to_string(), fastq.to_str().unwrap()])
        .args(&args)
        .stdout(stdout)
        .output()
        .expect("Failed to execute bwa");

    log_output(params, &output.stderr, "bwa mem complete; Running samtools sort");

    run_samtools_sort(params, align_sam, align_sorted_sam)?;
    run_samtools_view(params, align_sorted_sam, align_sorted_bam)
}

fn run_samtools_sort(params: &Paramsm, input: &Path, output: &Path) -> Result<(), Box<dyn Error>> {
    let command = "samtools";
    let args = &["sort", "-@", &params.threads.to_string(), "-o", output.to_str().unwrap(), input.to_str().unwrap()];

    let output = Command::new(command)
        .args(args)
        .output()
        .expect("Failed to execute samtools sort");

    log_output(params, &output.stderr, "Samtools sort complete; Running samtools view");

    Ok(())
}

fn run_samtools_view(params: &Paramsm, input: &Path, output: &Path) -> Result<(), Box<dyn Error>> {
    let command = "samtools";
    let args = &["view", "-b", "-@", &params.threads.to_string(), "-o", output.to_str().unwrap(), input.to_str().unwrap()];

    let output = Command::new(command)
        .args(args)
        .output()
        .expect("Failed to execute samtools view");

    log_output(params, &output.stderr, "Samtools view complete; Running samtools index");

    Ok(())
}

fn index_bam(params: &Paramsm, align_sorted_bam: &Path) -> Result<(), Box<dyn Error>> {
    let output = Command::new("samtools")
        .args(&["index", "-@", &params.threads.to_string(), align_sorted_bam.to_str().unwrap()])
        .output()
        .expect("Failed to execute samtools index");

    log_output(params, &output.stderr, "");

    Ok(())
}

fn cleanup(params: &Paramsm, align_sam: &Path, align_sorted_sam: &Path, fastq: &Path) -> Result<(), Box<dyn Error>> {
    fs::remove_file(align_sam)?;
    fs::remove_file(align_sorted_sam)?;
    fs::remove_file(fastq)?;

    if params.aligner.flavor == AlignerFlavor::STAR {
        let possible_tmp = params.output_path.join("_STARtmp");
        if possible_tmp.exists() {
            fs::remove_dir_all(possible_tmp)?;
        }
    }

    Ok(())
}

fn log_info(params: &Paramsm, message: &str) {
    if params.verbose {
        eprintln!("{}", message);
    }
    info!("{}", message);
}

fn log_output(params: &Paramsm, stderr: &[u8], message: &str) {
    if params.verbose {
        eprintln!("{}", String::from_utf8_lossy(stderr));
        eprintln!("{}", message);
    }
    info!("{}", String::from_utf8_lossy(stderr));
    info!("{}", message);
}

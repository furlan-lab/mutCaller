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


#[cfg(not(feature = "paris"))]
use log::*;

#[derive(Deserialize, Clone, PartialEq)]
pub struct Aligner {
    pub flavor: AlignerFlavor,
    pub loc: String,
    pub args: Option<Vec<String>>,
}

#[derive(Deserialize, Clone, PartialEq)]
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

pub fn align(params: &Paramsm) -> Result<(), Box<dyn Error>> {
    let align_sam = params.output_path.join("Aligned.sam");
    let align_sorted_sam = params.output_path.join("Aligned.sorted.sam");
    let align_sorted_bam = params.output_path.join("Aligned.sortedByCoord.out.bam");
    let fastq = params.output_path.join("mutcaller_R1.fq.gz");
    let outfolder = params.output_path.join("");
    let zcat_cmd = if cfg!(target_os = "macos") { "zcat <" } else { "zcat" };

    match params.aligner.flavor {
        AlignerFlavor::Minimap2 => {
            let args = params.aligner.args.clone();
            if params.verbose {
                eprintln!("Aligning reads using minimap2");
            }
            info!("Aligning reads using minimap2");
            let output = Command::new(&params.aligner.loc)
                .args(["--MD", "-Y", "-a", &params.genome, "-t", &params.threads.to_string(), "-o", align_sam.to_str().unwrap()])
                .args(args.unwrap_or_default())
                .arg(fastq.to_str().unwrap())
                .stderr(Stdio::piped())
                .stdout(Stdio::piped())
                .output()
                .expect("\n\n*******Failed to execute minimap2*******\n\n");
            if params.verbose {
                eprintln!("{}", String::from_utf8_lossy(&output.stderr));
                eprintln!("Minimap2 complete; Running samtools sort");
            }
            info!("{}", String::from_utf8_lossy(&output.stderr));
            info!("Minimap2 complete; Running samtools sort");
        }
        AlignerFlavor::STAR => {
            let args = params.aligner.args.clone();
            if params.verbose {
                eprintln!("Aligning reads using STAR");
            }
            info!("Aligning reads using STAR");
            let output = Command::new(&params.aligner.loc)
                .arg("--outFileNamePrefix").arg(outfolder.to_str().unwrap())
                .arg("--genomeDir").arg(&params.genome)
                .arg("--readFilesIn").arg(fastq.to_str().unwrap())
                .arg("--readNameSeparator").arg("space")
                .arg("--runThreadN").arg(params.threads.to_string())
                .arg("--outSAMunmapped").arg("Within").arg("KeepPairs")
                .arg("--outSAMtype").arg("BAM").arg("SortedByCoordinate")
                .arg("--outSAMattributes").arg("All")
                .arg("--readFilesCommand").arg(zcat_cmd)
                .args(args.unwrap_or_default())
                .stderr(Stdio::piped())
                .stdout(Stdio::piped())
                .output()
                .expect("\n\n*******Failed to execute STAR*******\n\n");
            if params.verbose {
                eprintln!("{}", String::from_utf8_lossy(&output.stderr));
                eprintln!("STAR complete; running samtools index");
            }
        }
        AlignerFlavor::BWA => {
            let stdout = File::create(align_sam.to_str().unwrap()).expect("failed to open output file for bwa mem");
            let args = params.aligner.args.clone();
            if params.verbose {
                eprintln!("Aligning reads using bwa mem");
            }
            info!("Aligning reads using bwa mem");
            let output = Command::new(&params.aligner.loc)
                .arg("mem")
                .arg("-t").arg(params.threads.to_string())
                .args(args.unwrap_or_default())
                .arg(&params.genome)
                .arg(fastq.to_str().unwrap())
                .stderr(Stdio::piped())
                .stdout(stdout)
                .output()
                .expect("\n\n*******Failed to execute bwa*******\n\n");
            if params.verbose {
                eprintln!("{}", String::from_utf8_lossy(&output.stderr));
                eprintln!("bwa mem complete; Running samtools sort");
            }
            info!("{}", String::from_utf8_lossy(&output.stderr));
            info!("bwa mem complete; Running samtools sort");
        }
    }

    let output = Command::new("samtools")
        .arg("sort")
        .arg("-@").arg(params.threads.to_string())
        .arg("-o").arg(align_sorted_sam.to_str().unwrap())
        .arg(align_sam.to_str().unwrap())
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .output()
        .expect("\n\n*******Failed to execute samtools view*******\n\n");
    if params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
        eprintln!("Samtools sort complete; Running samtools view");
    }
    info!("{}", String::from_utf8_lossy(&output.stderr));
    info!("Samtools sort complete; Running samtools view");

    let output = Command::new("samtools")
        .arg("view")
        .arg("-b")
        .arg("-@").arg(params.threads.to_string())
        .arg("-o").arg(align_sorted_bam.to_str().unwrap())
        .arg(align_sorted_sam.to_str().unwrap())
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .output()
        .expect("\n\n*******Failed to execute samtools sort*******\n\n");
    if params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
        eprintln!("Samtools view complete; Running samtools index");
    }
    info!("{}", String::from_utf8_lossy(&output.stderr));
    info!("Samtools view complete; Running samtools index");

    let output = Command::new("samtools")
        .arg("index")
        .arg("-@").arg(params.threads.to_string())
        .arg(align_sorted_bam.to_str().unwrap())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("\n\n*******Failed to execute samtools index*******\n\n");
    if params.verbose {
        eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    }
    info!("{}", String::from_utf8_lossy(&output.stderr));

    if !params.keep {
        match params.aligner.flavor {
            AlignerFlavor::Minimap2 => {
                fs::remove_file(align_sorted_sam.to_str().unwrap())?;
                fs::remove_file(align_sam.to_str().unwrap())?;
                fs::remove_file(fastq.to_str().unwrap())?;
            }
            AlignerFlavor::STAR => {
                fs::remove_file(fastq.to_str().unwrap())?;
                let possible_tmp = params.output_path.join("_STARtmp");
                if possible_tmp.exists() {
                    fs::remove_dir_all(possible_tmp)?;
                }
            }
            _ => (),
        }
    }

    Ok(())
}

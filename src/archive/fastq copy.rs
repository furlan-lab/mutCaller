/*

time ~/develop/mutCaller/target/release/fastq1 -t 1 --barcodes_file /Users/sfurlan/develop/mutCaller/data/737K-august-2016.txt.gz --fastq1 sequencer_R1.fastq.gz --fastq2 sequencer_R2.fastq.gz | gzip > out1.fq.gz
real    0m1.219s
echo $(zcat < tests/out1.fq.gz | wc -l)/4|bc

#compare to original mutcaller
cd ~/develop/mutCaller
time ~/develop/mutCaller/mutcaller -u -U 10 -b ~/develop/mutCaller/tests/sequencer_R1.fastq.gz -t ~/develop/mutCaller/tests/sequencer_R2.fastq.gz -l ~/develop/mutCaller/data/737K-august-2016.txt.gz &&
gzip tests/fastq_processed/sample_filtered.fastq



## full run of pipeline on testdata STAR
ml R/4.1.1-foss-2020b
ml SAMtools
export transcriptome=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A
cd ~/develop/mutCaller/tests
time ~/develop/mutCaller/target/release/fastq1 -t 1 --fastq1 sequencer_R1.fastq.gz --fastq2 sequencer_R2.fastq.gz --barcodes_file /home/sfurlan/develop/mutCaller/data/737K-august-2016.txt.gz | gzip > out1.fq.gz
zcat out1.fq.gz | head
/app/software/CellRanger/7.0.1/lib/bin/STAR --genomeDir /fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A/star --readFilesIn <(gunzip -c out1.fq.gz) \
  --runThreadN 1 --outSAMunmapped Within KeepPairs --outSAMtype BAM SortedByCoordinate
samtools view Aligned.sortedByCoord.out.bam | head
Rscript ~/develop/mutCaller/scripts/quantReads.R
time ~/develop/mutCaller/target/release/addtag -j _ --ibam Aligned.sortedByCoord.out.bam --obam Aligned.sortedByCoord.out.tagged.bam
samtools view Aligned.sortedByCoord.out.tagged.bam | head
time ~/develop/mutCaller/target/release/count -s _ -t 24 --ibam=Aligned.sortedByCoord.out.tagged.bam > counts_s.txt
sort -n -k3 -k2 -k1 counts_s.txt | uniq -c | sort -k2 -k3 -k4 > counts.sorted_s.txt

# real  0m23.133s

## full run of pipeline on testdata kallisto
ml kallisto/0.48.0-foss-2020b
export kallisto_index=/fh/scratch/delete90/furlan_s/targ_reseq/220819/kallisto_hg38_MYBindel
kallisto quant -i $kallisto_index -o kquant -l 400 -s 50 -t 18 --single --pseudobam out1.fq.gz
samtools view kquant/pseudoalignments.bam | head -n 1000
time ~/develop/mutCaller/target/release/count -t 24 --ibam=kquant/pseudoalignments.bam -a kallisto > counts_k.txt
sort -n -k3 -k2 -k1 counts_k.txt | uniq -c | sort -k2 -k3 -k4 > counts.sorted_k.txt
gzip counts.sorted_k.txt

## full run of pipeline on testdata using mm2
ml minimap2/2.24-GCCcore-11.2.0
ml SAMtools/1.16.1-GCC-11.2.0
cd ~/develop/mutCaller/tests
mkdir mm2
export fa=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
time ~/develop/mutCaller/target/release/fastq -t 1 --fastq1 sequencer_R1.fastq.gz --fastq2 sequencer_R2.fastq.gz --barcodes_file /home/sfurlan/develop/mutCaller/data/737K-august-2016.txt.gz | gzip > out1.fq.gz
zcat out1.fq.gz | head
minimap2 -a $fa --MD out1.fq.gz | samtools view -b -o mm2/Aligned.out.bam
samtools sort -o mm2/Aligned.out.sorted.bam mm2/Aligned.out.bam
samtools view mm2/Aligned.out.sorted.bam | head
samtools index mm2/Aligned.out.sorted.bam
time ~/develop/mutCaller/target/release/count -t 24 --ibam=mm2/Aligned.out.sorted.bam -v variants.tsv > counts_mm.txt
sort -n -k4 -k3 -k2 -k1 counts_mm.txt | uniq -c | sort -k2 -k3 -k4 > counts.sorted_mm.txt
gzip counts.sorted_mm.txt

*/





#[macro_use]
extern crate simple_log;
extern crate clap;
extern crate fastq;
extern crate flate2;

use simple_log::LogConfigBuilder;
use fastq::{parse_path, Record, RefRecord, each_zipped};
use clap::{App, load_yaml};
use flate2::{read};
use std::{
    // error::Error,
    ffi::OsStr,
    fs::File,
    io::{self, BufReader, BufRead},
    path::Path
};


// #[derive(Debug)]

fn lines_from_file(filename: &str) -> Vec<String> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(_why) => panic!("couldn't open {}", path.display()),
        Ok(file) => file,
    };
    if path.extension() == Some(OsStr::new("gz")){
        let buf = BufReader::new(read::GzDecoder::new(file));
        buf.lines()
            .map(|l| l.expect("Could not parse line"))
            .collect()
    }else{
        let buf = BufReader::new(file);
        buf.lines()
            .map(|l| l.expect("Could not parse line"))
            .collect()
    }
}


struct Params {
    fastq1: String,
    fastq2: String,
    ofastq: String,
    bcs: String,
    umi_len: u8,
    cb_len: u8,
    threads: usize,
    // max_reads: usize,
    // stop: bool,
    // debug: bool,
    name_sep: String,
}

fn load_params() -> Params {
    // let max_r = 0usize;
    // let stop = true;
    let yaml = load_yaml!("params_fastq1.yml");
    let params = App::from_yaml(yaml).get_matches();
    // let debug = params.value_of("debug").unwrap_or("false");
    // let debug = debug.to_string().parse::<bool>().unwrap();
    let fastq1 = params.value_of("fastq1").unwrap();
    let fastq2 = params.value_of("fastq2").unwrap();
    let ofastq = params.value_of("outfastq").unwrap_or("out.fastq.gz");
    let bcs = params.value_of("barcodes_file").unwrap_or("/Users/sfurlan/develop/mutCaller/data/737K-august-2016.txt.gz");
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<u8>().unwrap();
    let umi_len = params.value_of("umi_len").unwrap_or("10");
    let umi_len = umi_len.to_string().parse::<u8>().unwrap();
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<u8>().unwrap();
    let name_sep = params.value_of("name_sep").unwrap_or("|BARCODE=");
    Params{
        fastq1: fastq1.to_string(),
        fastq2: fastq2.to_string(),
        ofastq: ofastq.to_string(),
        bcs: bcs.to_string(),
        threads: threads as usize,
        umi_len: umi_len as u8,
        cb_len: cb_len as u8,
        name_sep: name_sep.to_string(),
    }
}

fn main() {
    let config = LogConfigBuilder::builder()
        .path("./fastq1.log")
        .size(1 * 100)
        .roll_count(10)
        .time_format("%Y-%m-%d %H:%M:%S.%f") //E.g:%H:%M:%S.%f
        .level("debug")
        .output_file()
        .build();

    let _ = simple_log::new(config);
    info!("starting!");
    let params = load_params();
    eprintln!("Running with {} thread(s)!", &params.threads);
    fastq(&params);
    info!("done!");
}


fn remove_whitespace(s: &mut String) {
    s.retain(|c| !c.is_whitespace());
}

// fn remove_whitespace_str(s: &str) -> String {
//     s.chars().filter(|c| !c.is_whitespace()).collect()
// }

fn fastq(params: &Params) {
    let mut cbvec = lines_from_file(&params.bcs);
    cbvec.sort_unstable();
    let _zip = true;
    let mut total_count: usize = 0;
    let mut nfound_count: usize = 0;
    let mut mmcb_count: usize = 0;
    let split_at = &params.umi_len + &params.cb_len;
    // let sep: Vec::<u8> = params.name_sep.as_bytes().to_vec();

    let fastq1 = &params.fastq1;
    let fastq2 = &params.fastq2;
    let _counts = (0u64, 0u64);
    let path = Path::new(&params.ofastq);
    let _file = match File::create(&path) {
        Err(_why) => panic!("couldn't open {}", path.display()),
        Ok(file) => file,
    };
    let mut writer = io::stdout();
    parse_path(Some(fastq1), |parser1| {
        parse_path(Some(fastq2), |parser2| {
            each_zipped(parser1, parser2, |rec1, rec2| {
                if rec1.is_some() & rec2.is_some(){
                    let r1 = &rec1.unwrap();
                    let r2 = &rec2.unwrap();
                    if r1.seq().contains(&b"N"[0]) | r2.seq().contains(&b"N"[0]){
                        nfound_count += 1;
                        total_count +=1;
                    }else{
                        total_count +=1;
                        let (barcode, _seq) = &r1.seq().split_at(split_at.into());
                        let (cb, _seq) = barcode.split_at(params.cb_len as usize);
                        match cbvec.binary_search(&std::str::from_utf8(cb).unwrap().to_string()) {
                            Ok(_u) => {
                                let mut readout = RefRecord::to_owned_record(&r2);
                                let _some_x = vec![b" "];
                                let mut new_header = std::str::from_utf8(&readout.head()).unwrap().to_string();
                                remove_whitespace(&mut new_header);
                                let _ = new_header.push_str(&params.name_sep);
                                let _ = new_header.push_str(&std::str::from_utf8(&barcode).unwrap().to_string());
                                readout.head = new_header.as_bytes().to_vec();

                                let _ = readout.write(&mut writer);
                            }
                            Err(_e) => {
                                mmcb_count +=1;
                            }
                        }
                    }
                }
                (true, true)
            })
            .expect("Invalid record.");
        })
        .expect("Unknown format for file 2.");
    })
    .expect("Unknown format for file 1.");
    eprintln!("Total number of reads processed: {}, {} of these had Ns, {} of these had BC not in whitelist", total_count, nfound_count, mmcb_count);
}

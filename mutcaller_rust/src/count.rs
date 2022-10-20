/**

cd ~/develop/mutCaller/mutcaller_rust/tests
time ~/develop/mutCaller/mutcaller_rust/target/release/count -t 24 --ibam=kquant/pseudoalignments.bam -a kallisto > counts_k.txt


time ~/develop/mutCaller/mutcaller_rust/target/release/count -t 24 --ibam=Aligned.sortedByCoord.out.tagged.bam > counts_s.txt


time ~/develop/mutCaller/mutcaller_rust/target/release/count -t 24 --ibam=Aligned.sortedByCoord.out.tagged.bam > counts_s.txt


##varianst
~/develop/mutCaller/mutcaller_rust/target/release/count -t 24 --ibam=Aligned.sortedByCoord.out.tagged.bam
**/

extern crate clap;
extern crate bam;

use std::io;
use clap::{App, load_yaml};
use std::str;
use crate::bam::RecordWriter;





// #[derive(Clone)]

struct Params {
    ibam: String,
    aligner: String,
    variants: String,
    split: String,
    joiner: String,
    threads: usize,
    cb_len: usize, 
    // umi_len: usize,
    read_len: usize,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let ibam = params.value_of("input_bam").unwrap();
    let aligner = params.value_of("aligner").unwrap_or("STAR");
    eprintln!("opening: {}", ibam.to_string());
    let variants = params.value_of("variants").unwrap_or("6,135195908,135195908");
    let joiner = params.value_of("joiner").unwrap_or(":");
    let split = params.value_of("split").unwrap_or("|BARCODE=");
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<usize>().unwrap() - 1;
    let read_len = params.value_of("read_len").unwrap_or("90");
    let read_len = read_len.to_string().parse::<usize>().unwrap();
    // let umi_len = params.value_of("umi_len").unwrap_or("10");
    // let umi_len = umi_len.to_string().parse::<usize>().unwrap() - 1;
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap() - 1;
    Params{
        ibam: ibam.to_string(),
        aligner: aligner.to_string(),
        variants: variants.to_string(),
        threads: threads,
        split: split.to_string(),
        joiner: joiner.to_string(),
        cb_len: cb_len,
        read_len: read_len,
    }
}

fn main() {
    let params = load_params();
    if params.variants=="6,135195908,135195908" {
        count_variants(&params)
    }
    if params.aligner == "kallisto"{
        count_kallisto(&params);
        return;
    }
    if params.aligner == "STAR"{
        count_star(&params);
        return;
    }
    
}


fn count_star(params: &Params) {
    let mut total: usize = 0;
    let mut goodreadcount: usize = 0;
    let (_read_threads, _write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as u16, ((*&params.threads/2) -1) as u16)
    } else {
        (0 as u16, 0 as u16)
    };
    let reader = bam::BamReader::from_path(params.ibam.to_string(), 0).unwrap();
    let _output = io::BufWriter::new(io::stdout());
    let header = reader.header().clone();
    let data = header.reference_names();
    let mut seqnames = Vec::new();
    for seq in data {
        seqnames.push(seq)
    }
    for record in reader {
        total += 1;
        let newrecord = record.unwrap();
        let seqname = match str::from_utf8(&newrecord.name()) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
        let cbumi= seqname.split(&params.split).nth(1).unwrap().to_string();
        let _modified_name = seqname.replace(&params.split, &params.joiner);
        let (cb_umi_s1, cb_umi_s2) = cbumi.split_at(params.cb_len+1);
        let mut good_read = false;
        let cigarmatch = format!("{}M", *&params.read_len);
        let cigar = newrecord.cigar().to_string();
        if cigar == cigarmatch{
            good_read = true
        }
        if newrecord.mapq() < 255 as u8 {
            good_read = false
        }
        if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
            goodreadcount += 1;
            println!("{} {} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string(), newrecord.start());
        }
    }
    eprintln!("Completed; {} total reads processed!", &total);
    eprintln!("{} good reads counted!", &goodreadcount);
}



fn count_kallisto(params: &Params) {
    let mut total: usize = 0;
    let mut goodreadcount: usize = 0;
    let (_read_threads, _write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as u16, ((*&params.threads/2) -1) as u16)
    } else {
        (0 as u16, 0 as u16)
    };
    let reader = bam::BamReader::from_path(params.ibam.to_string(), 0).unwrap();
    let _output = io::BufWriter::new(io::stdout());
    let header = reader.header().clone();
    let data = header.reference_names();
    let mut seqnames = Vec::new();
    for seq in data {
        seqnames.push(seq)
    }
    for record in reader {
        total += 1;
        let newrecord = record.unwrap();
        let seqname = match str::from_utf8(&newrecord.name()) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
        let cbumi= seqname.split(&params.split).nth(1).unwrap().to_string();
        let _modified_name = seqname.replace(&params.split, &params.joiner);
        let (cb_umi_s1, cb_umi_s2) = cbumi.split_at(params.cb_len+1);
        let mut good_read = false;
        let cigarmatch = format!("{}M", *&params.read_len);
        let cigar = newrecord.cigar().to_string();
        if cigar == cigarmatch{
            good_read = true
        }
        if newrecord.mapq() < 255 as u8 {
            good_read = false
        }
        if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
            goodreadcount += 1;
            println!("{} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string());
        }
    }
    eprintln!("Completed; {} total alignments processed!", &total);
    eprintln!("{} good alignments counted!", &goodreadcount);
}

fn process_variant(params: &Params)->bam::Region{
    let region = bam::Region::new(6,135195908,135195908);
    return region;
}

// fn count_variants(params: &Params, bar: bar, mapq: maq_Q, baseq: baseq){
fn count_variants(params: &Params){
    //         // """
    //         // :param bam_fil: bam_file
    //         // :param bar: good barcodes dict
    //         // :param indel: SNP = 0 indel= int(window size) insertion or deletion
    //         // :param mapq: mapping quality
    //         // :param baseq: base quality
    //         // :return: two dicts 1) list of all barcodes at given var pos
    //         // 2) dict with ref and alt barcodes
    //         // """
    // let barcodes = defaultdict(list);
    // 
    let mut reader = bam::IndexedReader::from_path(&params.ibam).unwrap();
    let output = io::BufWriter::new(io::stdout());
    let mut writer = bam::SamWriter::build()
        .write_header(false)
        .from_stream(output, reader.header().clone()).unwrap();
    // let (bar, mapQ, baseq) = 0usize;
    let region = process_variant(&params);
    for record in reader.fetch(&region).unwrap() {
        let record = record.unwrap();
        writer.write(&record).unwrap();
    }

}
    // // let barUcodes = defaultdict(list);
    // // let bar_count = [("ref", vec![]), ("alt", vec![])].iter().cloned().collect::<HashMap<_,_>>();
    // for column in bam::Pileup::with_filter(&mut reader, |record| record.flag().no_bits(1796)) 
    // with!(pysam.AlignmentFile(bam_fil, "rb") as pile){
    //     // logger.info("processing variant: {} {}  {}  {}  {}  {}  {}".format(self.chrm, self.start, self.end, self.ref, self.alt, self.gene, self.event));
    //     for pileupcolumn in pile.pileup(self.chrm, (self.start - 1), self.start, true, "nofilter", 100000000) {
    //         for read in pileupcolumn.pileups {
    //             if read.alignment.has_tag("CB")&&read.alignment.has_tag("UB") {
    //                 if bar.iter().all(|&x| x != read.alignment.get_tag("CB")) {
    //                     continue;
    //                 }
    //                 if indel > 0 {
    //                         if read.is_refskip||i32::from(read.alignment.mapping_quality) < i32::from(mapq) {
    //                             continue;
    //                         }
    //                 } else {
    //                     if read.is_del||read.is_refskip||i32::from(read.alignment.query_qualities[read.query_position]) < i32::from(baseq)||i32::from(read.alignment.mapping_quality) < i32::from(mapq) {
    //                         continue;
    //                     }:"?"
    //                 }
    //                 let mut q_name = read.alignment.query_name;
    //                 let mut q_tag = read.alignment.get_tag("CB");
    //                 let mut qU_tag = read.alignment.get_tag("UB");
    //                 let mut qstring = ((q_tag + ":") + qU_tag);
    //                 barcodes["{}".format(q_name)].append(q_tag);
    //                 barUcodes["{}".format(q_name)].append(qstring);
    //                 if indel > 0 {
    //                     q_name = read.alignment.query_name;
    //                     q_tag = read.alignment.get_tag("CB");
    //                     qU_tag = read.alignment.get_tag("UB");
    //                     qstring = ((q_tag + ":") + qU_tag);
    //                     barcodes["{}".format(q_name)].append(q_tag);
    //                     barUcodes["{}".format(q_name)].append(qstring);
    //                     if self.alt == "-" {
    //                     if read.alignment.has_tag("CB")&&read.query_position == None {
    //                     let mut alt_tag = read.alignment.get_tag("CB");
    //                     let mut altU_tag = read.alignment.get_tag("UB");
    //                     let mut ustringa = ((alt_tag + ":") + altU_tag);
    //                     bar_count["alt"].append(ustringa);
    //                     } else {
    //                     let mut ref_tag = read.alignment.get_tag("CB");
    //                     let mut refU_tag = read.alignment.get_tag("UB");
    //                     let mut ustring = ((ref_tag + ":") + refU_tag);
    //                     bar_count["ref"].append(ustring);
    //                     }
    //                     } else {
    //                     if read.alignment.has_tag("CB")&&read.indel == indel&&read.alignment.query_alignment_sequence[(read.query_position + 1)..((read.query_position + read.indel) + 1)] == self.alt {
    //                     let mut alt_tag = read.alignment.get_tag("CB");
    //                     let mut altU_tag = read.alignment.get_tag("UB");
    //                     let mut ustringa = ((alt_tag + ":") + altU_tag);
    //                     bar_count["alt"].append(ustringa);
    //                     } else {
    //                     let mut ref_tag = read.alignment.get_tag("CB");
    //                     let mut refU_tag = read.alignment.get_tag("UB");
    //                     let mut ustring = ((ref_tag + ":") + refU_tag);
    //                     bar_count["ref"].append(ustring);
    //                     }
    //                     }
    //                 } else {
    //                     if read.alignment.query_sequence[read.query_position] != self.alt {
    //                         let mut ref_tag = read.alignment.get_tag("CB");
    //                         let mut refU_tag = read.alignment.get_tag("UB");
    //                         let mut ustring = ((ref_tag + ":") + refU_tag);
    //                         bar_count["ref"].append(ustring);
    //                     } else {
    //                     if read.alignment.query_sequence[read.query_position] == self.alt {
    //                         let mut alt_tag = read.alignment.get_tag("CB");
    //                         let mut altU_tag = read.alignment.get_tag("UB");
    //                         let mut ustringa = ((alt_tag + ":") + altU_tag);
    //                         bar_count["alt"].append(ustringa);
    //                         }      
    //                     }
    //                 }
    //                 }
    //             }
    //         }
    //     }
    //     return (barUcodes, bar_count);



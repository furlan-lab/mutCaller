/**

cd ~/develop/mutCaller/tests
#kallisto
time ~/develop/mutCaller/target/release/count -t 24 --ibam=kquant/pseudoalignments.bam -a kallisto > counts_k.txt

#star
time ~/develop/mutCaller/target/release/count -t 24 --ibam=Aligned.sortedByCoord.out.tagged.bam > counts_s.txt

#mm2
time ~/develop/mutCaller/target/release/count -t 24 --ibam=mm2/Aligned.out.sorted.tagged.bam -v variants.tsv -m tags > counts.tsv
time ~/develop/mutCaller/target/release/count -t 24 --ibam=mm2/Aligned.out.sorted.tagged.bam -v variants.tsv -m header -s _ > counts.tsv

cat > variants.csv << EOL
seq,start,ref_nt,query_nt,name
chr12,112450407,A,G,PTPN11_227A>G
chr2,208248389,G,A,IDH1_132G>A
chr17,7674220,C,T,TP53_248C>T
chr6,135181325,A,G,MYBindel
EOL
cat variants.csv | sed 's/,/\t/g' > variants.tsv

chr6,135181308,135181308,gccagcaaggtgcatga,indel

samtools view mm2/Aligned.out.sorted.tagged.bam "chr6:135181308-135181308"
samtools view Aligned.sortedByCoord.out.tagged.bam "chr6:135181308-135181308"

samtools view mm2/Aligned.out.sorted.tagged.bam "chr12:112450407-112450407"




###TODO####
check sequencing quality of base in question
parameterize mapQ for mm2



~/develop/mutCaller/target/release/count -t 24 --ibam=mm2/Aligned.out.sorted.bam -v variants.tsv

**/
extern crate csv;
extern crate clap;
extern crate bam;
extern crate serde;
extern crate itertools;

// use std::fmt;
use std::io;
use clap::{App, load_yaml};
use std::str;
use std::error::Error;
use serde::Deserialize;
use std::fmt; 
use csv::ReaderBuilder;
// use csv::Reader;
use std::fs::File;
use std::io::{BufRead, BufReader};
use itertools::Itertools;
use flate2::GzBuilder;
use flate2::Compression;
use std::io::Write;


// #[derive(Clone)]

#[derive(Deserialize)]
struct Variant {
    seq: String,
    start: String,
    ref_nt: char,
    query_nt: char,
    name: String,
}

// Implement `Display` for `Variant`.
impl fmt::Display for Variant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use `self.number` to refer to each positional data point.
        write!(f, "seq: {} start: {} ref_nt: {} query_nt: {} name: {}", self.seq, self.start, self.ref_nt, self.query_nt, self.name)
    }
}


// Implement `Display` for `Record`.
impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use `self.number` to refer to each positional data point.
        write!(f, "cb: {} umi: {} seq: {} loc: {} variant: {} call: {}", self.cb, self.umi, self.seq, self.loc, self.variant, self.call)
    }
}



struct Params {
    ibam: String,
    aligner: String,
    variants: String,
    split: String,
    joiner: String,
    threads: usize,
    cb_len: usize, 
    read_len: usize,
    method: String,
}



fn read_csv(params: &Params) -> Result<Vec<Variant>, Box<dyn Error>> {
    // Build the CSV reader and iterate over each record.
//     let data = "\
// seq\tstart\tref_nt\tquery_nt\tname
// chr12\t112450407\tA\tG\tPTPN11_227A>G
// chr2\t208248389\tG\tA\tIDH1_132G>A
// chr17\t7674220\tC\tT\tTP53_248C>T";
    eprintln!("opening variants file: {}", &params.variants.to_string());
    let file = File::open(&params.variants.to_string()).unwrap();
    let reader = BufReader::new(file);
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(reader);
    let mut csvdata = Vec::new();
    for result in rdr.deserialize() {
        let record: Variant = result?;
        csvdata.push(record);
    }
    // let dummyvariant: Variant = Variant {seq: String::from("chr12"),
    //         start: String::from("112450407"),
    //         ref_nt: 'A',
    //         query_nt: 'G',
    //         name: String::from("PTPN11_227A>G"),
    //     };
    // csvdata.push(dummyvariant);
    Ok(csvdata)
}


fn load_params() -> Params {
    let yaml = load_yaml!("params_count.yml");
    let params = App::from_yaml(yaml).get_matches();
    let ibam = params.value_of("input_bam").unwrap();
    let aligner = params.value_of("aligner").unwrap_or("mm2");
    let variantstring = params.value_of("variants").unwrap();
    let joiner = params.value_of("joiner").unwrap_or(":");
    let split = params.value_of("split").unwrap_or("|BARCODE=");
    let cb_len = params.value_of("cb_len").unwrap_or("16");
    let cb_len = cb_len.to_string().parse::<usize>().unwrap() - 1;
    let read_len = params.value_of("read_len").unwrap_or("90");
    let read_len = read_len.to_string().parse::<usize>().unwrap();
    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap() - 1;
    let method = params.value_of("method").unwrap_or("header");
    Params{
        ibam: ibam.to_string(),
        aligner: aligner.to_string(),
        variants: variantstring.to_string(),
        threads: threads,
        split: split.to_string(),
        joiner: joiner.to_string(),
        cb_len: cb_len,
        read_len: read_len,
        method: method.to_string(),
    }
}

fn main() {
    let params = load_params();
    if params.aligner == "mm2" {
        let csvdata = read_csv(&params).unwrap();
        // let f = File::create("counts_mm.txt.gz")?;
        // let mut gz = GzBuilder::new()
        //                 .filename("counts_mm.txt")
        //                 .write(f, Compression::default());
        // for (count, record) in cdata {
        //     let count_str = record+&" ".to_owned()+&(count.to_string());
        //     gz.write_all(count_str.as_bytes())?;
        // }
        // gz.finish()?;
        let mut count_vec = Vec::new();
        for variant in csvdata {
            eprintln!("Processing variant: {}", variant);
            eprintln!("opening bam: {}", &params.ibam.to_string());
            if params.method.to_string() == "tags"{
                count_vec.push(count_variants_mm(&params, variant));
                // let cdata = count(&sdata);

            } else {
                count_vec.push(count_variants_mm2(&params, variant));
                // let cdata = count(&sdata);
            }
            
        }
        let _none = writer_fn(count_vec, "counts_mm.txt.gz".to_string());
        return;
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

fn writer_fn (count_vec: Vec<Vec<Vec<u8>>>, fname: String) -> Result<(), Box<dyn Error>> {
        let f = File::create(fname)?;
        let mut gz = GzBuilder::new()
                        .filename("counts_mm.txt.gz")
                        .write(f, Compression::default());
        for result in count_vec {
            for line in result {
                gz.write_all(&line)?;
            }
        }
        gz.finish()?;
        Ok(())
}



#[derive(Deserialize, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct Record {
    cb: String,
    umi: String,
    seq: String,
    loc: usize,
    variant: String,
    call: String,
}


fn count_star(params: &Params) {
    eprintln!("Counting star reads");
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
        if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
            goodreadcount += 1;
            println!("{} {} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string(), newrecord.start());
        }
    }
    eprintln!("Completed; {} total reads processed!", &total);
    eprintln!("{} good reads counted!", &goodreadcount);
}



fn count_kallisto(params: &Params) {
    eprintln!("Counting kallisto reads");
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
        if good_read && ((newrecord.flag().to_string()=="Flag(16)") | (newrecord.flag().to_string()=="Flag(0)")){
            goodreadcount += 1;
            println!("{} {} {}", cb_umi_s1, cb_umi_s2, seqnames[newrecord.ref_id() as usize].to_string());
        }
    }
    eprintln!("Completed; {} total alignments processed!", &total);
    eprintln!("{} good alignments counted!", &goodreadcount);
}


fn process_variant(ref_id: u32, start: u32)->bam::Region{
    let region = bam::Region::new(ref_id,start - 1,start - 1);
    return region;
}


fn count_variants_mm2(params: &Params, variant: Variant) -> Vec<Vec<u8>>{
    eprintln!("Processing using cb and umi in header");
    let mut total: usize = 0;
    let seqname = variant.seq;
    let start = variant.start.parse::<u32>().unwrap();
    let vname = variant.name;
    let mut reader = bam::IndexedReader::build()
        .additional_threads(*&params.threads as u16)
        .from_path(&params.ibam).unwrap();
    let mut seqnames = Vec::new();
    let mut result = "";
    let query_nt = variant.query_nt as char;
    let header = reader.header().clone();
    let hdata = header.reference_names();
    for seq in hdata {
        seqnames.push(seq)
    }
    let mut data = Vec::new();
    let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
    let region = process_variant(ref_id as u32, start);
    for record in reader.fetch_by(&&region, |record| record.mapq() >= 4 && (record.flag().all_bits(0 as u16) || record.flag().all_bits(16 as u16))).unwrap(){
        total+=1;
        let readheader = match str::from_utf8(record.as_ref().unwrap().name()) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
        let cbumi= readheader.split(&params.split).nth(1).unwrap().to_string();
        let (cb, umi) = cbumi.split_at(params.cb_len+1);
        for entry in record.as_ref().unwrap().alignment_entries().unwrap() {
            if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
                if region.start() == ref_pos {
                    if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
                        if ref_nt as char == record_nt as char {
                            result = "ref";
                        } else if record_nt as char == query_nt{
                            result = "query";
                        } else {
                            result = "other";
                        }
                            data.push(format!("{} {} {} {} {} {}", &cb, &umi, seqname, ref_pos, vname, result))
                        }
                    } else {
                        continue
                    }
            } else {
                continue
            }        }
    }
    eprintln!("Found {} reads spanning this variant!", total);
    data.sort();
    let mut out_vec = Vec::new();
    let cdata = data.into_iter().dedup_with_count();
    for (count, record) in cdata {
       let count_str = record+&" ".to_owned()+&(count.to_string()+&"\n".to_owned());
        out_vec.push(count_str.as_bytes().to_owned());
    }
    return out_vec;
}


fn count_variants_mm(params: &Params, variant: Variant) -> Vec<Vec<u8>>{
    eprintln!("Processing using cb and umi in BAM tags");
    let mut total: usize = 0;
    let seqname = variant.seq;
    let start = variant.start.parse::<u32>().unwrap();
    let vname = variant.name;
    let mut reader = bam::IndexedReader::build()
        .additional_threads(*&params.threads as u16)
        .from_path(&params.ibam).unwrap();
    let mut seqnames = Vec::new();
    let mut cb;
    let mut umi;
    let mut result = "null";
    let query_nt = variant.query_nt as char;
    let header = reader.header().clone();
    let hdata = header.reference_names();
    for seq in hdata {
        seqnames.push(seq)
    }
    let mut data = Vec::new();
    let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
    let region = process_variant(ref_id as u32, start);
    for record in reader.fetch_by(&&region, |record| record.mapq() > 4 && (record.flag().all_bits(0 as u16) || record.flag().all_bits(16 as u16))).unwrap(){
        total+=1;
        match record.as_ref().unwrap().tags().get(b"CB") {
            Some( bam::record::tags::TagValue::String(cba, _)) => {
                cb = str::from_utf8(&cba).unwrap().to_string();
            },
            _ => panic!("Unexpected type"),
        }
        match record.as_ref().unwrap().tags().get(b"UB") {
            Some( bam::record::tags::TagValue::String(uba, _)) => {
                umi = str::from_utf8(&uba).unwrap().to_string();
            },
            _ => panic!("Unexpected type"),
        }
        for entry in record.as_ref().unwrap().alignment_entries().unwrap() {
            if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
                if region.start() == ref_pos {

                    if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
                        if ref_nt as char == record_nt as char {
                            result = "ref";
                        } else if record_nt as char == query_nt{
                            result = "query";
                        } else {
                            result = "other";
                        }
                            data.push(format!("{} {} {} {} {} {}", &cb, &umi, seqname, ref_pos, vname, result))
                        }
                    } else {
                        continue
                    }
            } else {
                continue
            }        }
    }
    eprintln!("Found {} reads spanning this variant!", total);
    data.sort();
    let mut out_vec = Vec::new();
    let cdata = data.into_iter().dedup_with_count();
    for (count, record) in cdata {
        let count_str = record+&" ".to_owned()+&(count.to_string()+&"\n".to_owned());
        // let count_str = record+&" ".to_owned()+&(count.to_string())+&"\n".to_owned();
        out_vec.push(count_str.as_bytes().to_owned());
    }
    return out_vec;
}

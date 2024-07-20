use crate::mutcaller::Paramsm;
use crate::utils::lines_from_file;
use flate2::{GzBuilder, Compression};
use std::error::Error;
use std::fs::File;
use std::str;
use fastq::{parse_path, each_zipped, OwnedRecord, Record};
use simplelog::info;

pub fn replace_slice<T>(source: &mut [T], from: &[T], to: &[T])
where
    T: Clone + PartialEq,
{
    if source.starts_with(from) {
        source[..from.len()].clone_from_slice(to);
    }

    if source.len() > from.len() {
        replace_slice(&mut source[from.len()..], from, to);
    }
}

fn fix_fastq_header(mut header: Vec<u8>, split: &[u8], barcode: &[u8]) -> Vec<u8> {
    replace_slice(&mut header, &[32], &[95]);
    header.extend_from_slice(split);
    header.extend_from_slice(barcode);
    header
}

pub fn fastq(params: &Paramsm) -> Result<(), Box<dyn Error>> {
    let outfastq_temp = params.output_path.join("mutcaller_R1.fq.gz");
    let outfastq = outfastq_temp.to_str().unwrap();
    let split = "|BARCODE=".as_bytes();
    let mut cbvec = lines_from_file(&params.bcs);
    cbvec.sort_unstable();
    let mut total_count: usize = 0;
    let mut nfound_count: usize = 0;
    let mut mmcb_count: usize = 0;
    let mut mm_count: usize = 0;
    let split_at = params.umi_len + params.cb_len;

    let f = File::create(&outfastq)?;
    let mut writer = GzBuilder::new()
        .filename(outfastq)
        .write(f, Compression::default());
    parse_path(Some(&params.fastq1), |parser1| {
        parse_path(Some(&params.fastq2), |parser2| {
            each_zipped(parser1, parser2, |rec1, rec2| {
                if let (Some(r1), Some(r2)) = (rec1, rec2) {
                    if r1.seq().contains(&b'N') || r2.seq().contains(&b'N') {
                        nfound_count += 1;
                    } else if r2.seq().len() != r2.qual().len() {
                        mm_count += 1;
                    } else {
                        total_count += 1;
                        let (barcode, _) = r1.seq().split_at(split_at);
                        let (cb, _) = barcode.split_at(params.cb_len);
                        if cbvec.binary_search(&str::from_utf8(cb).unwrap().to_string()).is_ok() {
                            let new_header = fix_fastq_header(r2.head().to_vec(), split, barcode);
                            let readout = OwnedRecord {
                                head: new_header,
                                seq: r2.seq().to_vec(),
                                sep: Some(vec![43]),
                                qual: r2.qual().to_vec(),
                            };
                            readout.write(&mut writer).unwrap();
                        } else {
                            mmcb_count += 1;
                        }
                    }
                }
                (true, true)
            }).expect("Invalid record.");
        }).expect("Unknown format for file 2.");
    }).expect("Unknown format for file 1.");

    info!("\n\n\tTotal number of reads processed: {}, {} of these had Ns, {} of these had BC not in whitelist, {} of these had mismatched length of sequence and quality\n", total_count, nfound_count, mmcb_count, mm_count);
    if params.verbose {
        eprintln!("Total number of reads processed: {}, {} of these had Ns, {} of these had BC not in whitelist, {} of these had mismatched length of sequence and quality\n", total_count, nfound_count, mmcb_count, mm_count);
    }
    Ok(())
}

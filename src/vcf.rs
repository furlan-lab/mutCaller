use vcf::{VCFReader, U8Vec, VCFHeaderFilterAlt, VCFError, VCFRecord};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::BufReader;

fn read_vcf(file: String) -> Result<(), VCFError> {
    let mut reader = VCFReader::new(BufReader::new(MultiGzDecoder::new(File::open(
        file,
    )?)))?;

    // access FILTER contents
    assert_eq!(
        Some(VCFHeaderFilterAlt {
            id: b"PASS",
            description: b"All filters passed"
        }),
        reader.header().filter(b"PASS")
    );

    // access INFO contents
    assert_eq!(
        b"Stop position of the interval",
        reader.header().info(b"END").unwrap().description
    );

    // prepare VCFRecord object
    let mut vcf_record = reader.empty_record();

    // read one record
    reader.next_record(&mut vcf_record)?;

   // get record attributes
   assert_eq!(vcf_record.chromosome, b"13");
   assert_eq!(vcf_record.position, 32872836);
   assert_eq!(vcf_record.id, Vec::<U8Vec>::new());
   assert_eq!(vcf_record.reference, b"A");
   assert_eq!(vcf_record.alternative, vec![b"C"]);
   assert_eq!(vcf_record.qual, Some(495.23));
   assert_eq!(vcf_record.info(b"AC"), Some(&vec![b"1".to_vec()]));
   assert_eq!(
       vcf_record.genotype(b"SRP150637__HG00099", b"GT"),
       Some(&vec![b"0/0".to_vec()])
   );
   assert_eq!(
       vcf_record.genotype(b"SRP150637__HG00099", b"AD"),
       Some(&vec![b"31".to_vec(), b"0".to_vec()])
   );

    Ok(())
}

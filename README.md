<img width="150" alt="image" src="mutcaller.png">


#                       mutCaller
Pipeline for extracting variants from single cell genomics data


### Overview

mutCaller is a command line tool written in Rust for identifying SNVs in single cell data that contains a cell barcode (CB) and UMI such as 10X genomics data files.  With mutCaller, you can supply variants (VCF file support coming), and obtain a counts file of the number of variants that map to your query and their associated cell-barcodes and umis.  mutCaller using the takes as input either unaligned fastqs (UNALIGNED function) or an aligned BAM file with the CB and UMI as tags (ALIGNED function).

### Installation

mutCaller is written in Rust.  To install the rust compiler go to https://www.rust-lang.org/tools/install.  mutCaller requires two additional tools be available on the command line, minimap2 (https://github.com/lh3/minimap2) and samtools (https://samtools.github.io). 

To install mutCaller:
1. clone the repository by typing `https://github.com/furlan-lab/mutCaller.git` from the location you want to build from
2. enter the cloned repo by typing `cd mutCaller`
3. build by typing `cargo build --release`
4. the build process will create a self contained binary executable file in `targets/release` directory called `mutCaller`
5. move this binary elsewhere if desired (ideally somewhere referenced by your PATH environment variable - e.g. `~/.local/bin`)

### Updates

**version 0.30** - 5/9/23 - code cleanup and implemented bam option
**version 0.22** - 5/9/23 - first alpha release


### Usage

##### Variants File.

First users will generate a simple, headerless 'variants_file' that lists the variants to be counted (tsv.file). The alleles_file should look something like this:

```plaintext
seqname\tstart\tref_nt\tquery_nt\tname
```
**More detailed explanation:**
1. seqname - e.g. 'chr1', 'chr2', etc
2. position - 1-indexed position (e.g, '112450407')
3. ref_nt - nucleotide of the reference at this location
4. query_nt - nucleotide of your query
5. name - a string given to the name the variant in the outpull file

**three full lines should look something like this**

```plaintext
chr12   112450407   A   G   PTPN11_227A>G
chr12   208248389   G   A   IDH1_132G>A
chr17   7674220 C   T   TP53_248C>T
```

##### Barcodes file

A barcode whitelist.  v3 3prime and v2 5prime 10X whitelists are provided in the data folder.


##### Invocation

**To run mutCaller simply type:**

```plaintext
mutCaller ALIGNED --bam <file.bam> -a <variants.tsv> -o <folder>
mutCaller UNALIGNED --barcodes_file <barcodes_file> --fastq1 <fastq1> --fastq2 <fastq2> --genome <genome> --variants <variants.tsv>
```


##### Help menu

```plaintext
mutcaller 0.3.0
Scott Furlan
Single nucleotide variant counting pipeline for single cell genomics data

USAGE:
    mutcaller [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    ALIGNED      Count variants in previously aligned data
    UNALIGNED    Count variants after aligning data using minimap2
    help         Prints this message or the help of the given subcommand(s)
```

ALIGNED help
```plaintext

mutcaller-ALIGNED
Count variants in previously aligned data

USAGE:
    mutcaller ALIGNED [FLAGS] [OPTIONS] --bam <bam> --variants <variants>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information
    -q, --quiet      use this flag to run in quiet mode (no verbosity)

OPTIONS:
    -b, --bam <bam>              aligned bam file with cell barcode and umi in tags
    -c, --cb_tag <cb>            bam tag containing cell barcode; default = 'CB'
    -o, --output <output>        output filename for compressed counts (defaults to counts_mm.txt.gz)
    -t, --threads <threads>      threads
    -u, --umi_tag <umi>          bam tag containing umi; default = 'XM'
    -v, --variants <variants>    path to variants.tsv file (SNVs with mm2 only supported currently) with the following
                                 formatting per line - seqname\tstart\tref_nt\tquery_nt\tname; e.g.
                                 chr12,112450407,A,G,PTPN11_227A>G

```

UNALIGNED help

```plaintext
mutcaller-UNALIGNED
Count variants after aligning data using minimap2

USAGE:
    mutcaller UNALIGNED [FLAGS] [OPTIONS] --barcodes_file <barcodes_file> --fastq1 <fastq1> --fastq2 <fastq2> --genome <genome> --variants <variants>

FLAGS:
    -h, --help          Prints help information
    -k, --keep_files    use this flag to keep files (default is remove intermediate files)
    -V, --version       Prints version information
    -q, --quiet         use this flag to run in quiet mode (no verbosity)

OPTIONS:
    -a, --aligner <aligner>                aligner software - currently mm2 (default) and kallisto are supported
    -b, --barcodes_file <barcodes_file>    barcodes_file
    -c, --cb_length <cb_len>               length of umi sequence
    -i, --fastq1 <fastq1>                  input fastq with barcodes
    -j, --fastq2 <fastq2>                  input fastq with read
    -g, --genome <genome>                  fasta for minimap2 or transcriptome index for kallisto
    -o, --output <output>                  output filename for compressed counts (defaults to counts_mm.txt.gz)
    -r, --read_len <read_len>              read 2 length (default 90)
    -t, --threads <threads>                threads
    -u, --umi_length <umi_len>             length of umi sequence
    -v, --variants <variants>              path to variants.tsv file (SNVs with mm2 only supported currently) with the
                                           following formatting per line - seqname\tstart\tref_nt\tquery_nt\tname; e.g.
                                           chr12,112450407,A,G,PTPN11_227A>G

```









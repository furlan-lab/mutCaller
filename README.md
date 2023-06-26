<img width="300" alt="image" src="mutcaller.png">

Pipeline for extracting variants from single cell genomics data


### Overview

mutCaller is a command line tool written in Rust for identifying SNVs in single cell data that contains a cell barcode (CB) and UMI such as 10X genomics data files.  With mutCaller, you can supply variants (in a convenient tsv file or by using a VCF), and obtain a counts file of the number of variants that map to your query and their associated cell-barcodes and umis.  mutCaller using the takes as input either unaligned fastqs (UNALIGNED function) or an aligned BAM file with the CB and UMI as tags (ALIGNED function).

### Installation

mutCaller is written in Rust.  To install the rust compiler go to https://www.rust-lang.org/tools/install.  If using the UNALIGNED function mutCaller requires additional tools be available on the command line, minimap2 (https://github.com/lh3/minimap2) or STAR (https://github.com/alexdobin/STAR) and samtools (https://samtools.github.io). 

To install mutCaller:
1. clone the repository by typing `https://github.com/furlan-lab/mutCaller.git` from the location you want to build from
2. enter the cloned repo by typing `cd mutCaller`
3. build by typing `cargo build --release`
4. the build process will create a self contained binary executable file in `targets/release` directory called `mutCaller`
5. move this binary elsewhere if desired (ideally somewhere referenced by your PATH environment variable - e.g. `~/.local/bin`)

### Updates
- **version 0.5.\*** - 6/19/23 - indels supported, major code cleanup, checks variants.tsv file to ensure ref_nt matches the reference (for UNALIGNED-minimap2 runs only); added feature for passing additional arguments to aligner; fixed a bug that alllwed reads with different seq and qual lengths to tank STAR.
- **version 0.4.\*** - 6/14/23 - majorcode cleanup, added STAR aligner option, added vcf file support including a VARIANT subfunction that allows you to troubleshoot variants files
- **version 0.30** - 5/9/23 - code cleanup and implemented bam option
- **version 0.22** - 5/9/23 - first alpha release


### Usage

Note windows not supported; tested on MacOS and Linux (Ubuntu).

#### Variants file (user generated tsv)

Users can generate a simple, 'variants_file' that lists the variants to be counted (tsv.file). The variants file should look like this:

```plaintext
seqname\tstart\tref_nt\tquery_nt\tname
```
**More detailed explanation:**
1. seqname - e.g. 'chr1', 'chr2', etc
2. position - 1-indexed position (e.g, '112450407')
3. ref_nt - nucleotide(s) of the reference at this location
4. query_nt - nucleotide(s) of your query
5. name - a string given to the name the variant in the outpull file
6. for all indels, the first nt provided should be contained in the reference with the following characters representing the insertion or deletion

**three full lines should look something like this**

Note for duplications, it is hard to predict which the aligner will detect as an insertion.  We recommend for duplications that you include an insertion query assuming the first repeated sequence is detected as an insertion, then repeat for the second repeated sequence.  We have given an example below

```plaintext
seq     start  ref_nt query_nt name
chr12   112450407   A   G   PTPN11_227A>G
chr12   112450405   T   TGAG    PTPN11_insT>TGAG
chr12   112450407   AGTT    A   PTPN11_delAGTT>A
chr6    135195908   C   CGCCAGCAAGGTGCATG   MYBdup1
chr6    135195925   A   AGCCAGCAAGGTGCATG   MYBdup2
```

#### VCF support

Mutcaller will accept VCF files.  Not all versions have been tested.  With the VARIANTS function, a user can filter records from a VCF and optionally output a variants.tsv.  It is also possible to supply a VCF to the ALIGNED and UNALIGNED functions in mutcaller directly.  See examples [here](EXAMPLES.md)

#### Barcodes file

A barcode whitelist.  v3 3prime and v2 5prime 10X whitelists are provided in the data folder.


#### Invocation

**To run mutCaller simply type:**

```plaintext
mutCaller ALIGNED --bam <file.bam> -s <variants.tsv> -o <folder>
mutCaller UNALIGNED --barcodes_file <barcodes_file> --fastq1 <fastq1> --fastq2 <fastq2> --genome <genome> --variants <variants.tsv>
```

#### OTHER EXAMPLES
More examples on how to use mutcaller can be found [here](EXAMPLES.md)



### Help menus

```plaintext
mutcaller 0.5.4
Scott Furlan
Single nucleotide variant counting pipeline for single cell genomics data

USAGE:
    mutcaller [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    ALIGNED      Count variants in previously aligned data
    ARGPARSE     test arg parsing
    UNALIGNED    Count variants after aligning data
    VARIANTS     variants file debugging
    help         Prints this message or the help of the given subcommand(s)

Note: this is a work in progress. Use with caution.
```

##### ALIGNED help
```plaintext

mutcaller-ALIGNED
Count variants in previously aligned data

USAGE:
    mutcaller ALIGNED [FLAGS] [OPTIONS] --bam <bam> --variants <variants>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information
    -v, --verbose    use this flag to run in verbose mode

OPTIONS:
    -b, --bam <bam>              aligned, sorted bam file with cell barcode and umi in tags; see cb_tag/umi_tag
    -c, --cb_tag <cb>            bam tag containing cell barcode; default = 'CB'
    -o, --output <output>        output path (defaults to 'out'); inside folder a counts file will be called
                                 'counts.txt.gz', log will be called mutcaller.log
    -q, --qual <qual>            filter for variant quality (float); default = 95.0; only used if VCF file type is
                                 supplied
    -t, --threads <threads>      threads
    -u, --umi_tag <umi>          bam tag containing umi; default = 'XM'
    -s, --variants <variants>    path to variants.tsv or vcf file; For tsv, example formating =
                                 seqname\tstart\tref_nt\tquery_nt\tname; e.g. chr12,112450407,A,G,PTPN11_227A>G

```

##### UNALIGNED help

```plaintext

mutcaller-UNALIGNED
Count variants after aligning data

USAGE:
    mutcaller UNALIGNED [FLAGS] [OPTIONS] --barcodes_file <barcodes_file> --fastq1 <fastq1> --fastq2 <fastq2> --genome <genome> --variants <variants>

FLAGS:
    -h, --help          Prints help information
    -k, --keep_files    use this flag to keep files (default is remove intermediate files)
    -V, --version       Prints version information
    -v, --verbose       use this flag to run in verbose mode

OPTIONS:
    -x, --add_aligner_args <add_aligner_args>
            additional args supplied to aligner; MUST wrap in quotes; MUST also be supplied using the long flag with
            double dash and equal sign e.g. --aligner_args="--some_complex_arg parameter --some_other_arg
            other_parameter", i.e. CANNOT use short flag "-x" or long flag with a space e.g. "--aligner_args
            --some_complex_arg_parameter"
    -a, --aligner <aligner>
            aligner software - currently "minimap" (default), and "STAR" are supported; if not available on command
            line, supply in loc_aligner argument
    -l, --aligner_loc <aligner_loc>              path to aligner e.g. /app/software/CellRanger/6.0.1/lib/bin/STAR
    -b, --barcodes_file <barcodes_file>          barcodes_file
    -c, --cb_length <cb_len>                     length of umi sequence
    -i, --fastq1 <fastq1>                        input fastq with barcodes
    -j, --fastq2 <fastq2>                        input fastq with read
    -g, --genome <genome>
            fasta for minimap2 and genome index location for STAR; for checking variants, fasta must be indexed

    -o, --output <output>
            output path (defaults to 'out'); inside folder a counts file will be called 'counts.txt.gz', log will be
            called mutcaller.log
    -q, --qual <qual>
            filter for variant quality (float); default = 95.0; only used if VCF file type is supplied

    -r, --read_len <read_len>                    read 2 length (default 90)
    -t, --threads <threads>                      threads
    -u, --umi_length <umi_len>                   length of umi sequence
    -s, --variants <variants>
            path to variants.tsv or vcf file; For tsv, example formating = seqname\tstart\tref_nt\tquery_nt\tname; e.g.
            chr12,112450407,A,G,PTPN11_227A>G

```









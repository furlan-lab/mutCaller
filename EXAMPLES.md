<img width="200" alt="image" src="mutcaller.png">


#                       mutCaller working examples

### Build mutCaller and install binary

```sh
loc=~/develop/mutCaller # or location where you have cloned the repository
cd $loc  
cargo build --release
cp target/release/mutcaller ~/.local/bin
```

### Run UNALIGNED on short read fastqs using mm2

```sh
loc=~/develop/mutCaller # or location where you have cloned the repository
bc=$loc/data/737K-august-2016.txt.gz  #barcode whitelist
fa=/Users/sfurlan/refs/genome.fa #genome location i.e. GRCh38
#fa=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A/fasta/genome.fa
mutcaller UNALIGNED -t 8 -g $fa -b $bc -v $loc/tests/variants.tsv \
          --fastq1 $loc/tests/sequencer_R1.fastq.gz \
          --fastq2 $loc/tests/sequencer_R2.fastq.gz
```

### Run UNALIGNED on short read fastqs using STAR

```sh
loc=~/develop/mutCaller # or location where you have cloned the repository
bc=$loc/data/737K-august-2016.txt.gz  #barcode whitelist
fa=/fh/fast/furlan_s/grp/refs/GRCh38/refdata-gex-GRCh38-2020-A/star
../target/release/mutcaller UNALIGNED \
                        -t 8 -g $fa -b $bc -v variants.tsv -a STAR -l /app/software/CellRanger/6.0.1/lib/bin/STAR \
                        -o out_star --fastq1 sequencer_R1.fastq.gz \
                        --fastq2 sequencer_R2.fastq.gz



cd /fh/scratch/delete90/furlan_s/targ_reseq/230117_Sami/AML_1101_merge/align_nodedup
ml SAMtools/1.11-GCC-10.2.0
samtools view -bs 42.01 aln_sncr_fc_Malig.mdtag.bam > subsampleMalig.mdtag.bam

cd /fh/scratch/delete90/furlan_s/targ_reseq/230117_Sami/AML_1101_merge/align_nodedup

~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v ~/develop/mutCaller/tests/variants.tsv -t 1 -o test1
~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v ~/develop/mutCaller/tests/variants.tsv -t 4


~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 1
~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 22

~/develop/mutCaller/target/release/mutcaller ALIGNED -t 20 -b aln_sncr_fc_Malig.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv



~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 1
~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 22

sbatch -n 1 -c 1 -p campus-new -M gizmo --mem-per-cpu=21000MB --wrap='~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 1 -o out1'
sbatch -n 1 -c 30 -p campus-new -M gizmo --mem-per-cpu=21000MB --wrap='~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 30 -o out30'
sbatch -n 1 -c 1 -p campus-new -M gizmo --mem-per-cpu=21000MB --wrap='~/develop/mutCaller/target/release/mutcaller ALIGNED -b aln_sncr_fc_Malig.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 1 -o out1_real'


head /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv


cat > breaksit.csv << EOL
seq,start,ref_nt,query_nt,name
chr14,31113320,T,C,HECTD1_25831T>C
EOL


sed -E 's/("([^"]*)")?,/\2\t/g' breaksit.csv > breaksit.tsv


~/develop/mutCaller/target/release/mutcaller ALIGNED -b aln_sncr_fc_Malig.mdtag.bam -v breaksit.tsv -t 1 -o out1_break

samtools view aln_sncr_fc_Malig.mdtag.bam | grep 'm64272e_230122_062448/100729383/ccs/9569_11265' > breaksit.bam
samtools view -b aln_sncr_fc_Malig.mdtag.bam "chr14:31113320-31113321" > breaksit2.bam


~/develop/mutCaller/target/release/mutcaller ALIGNED -b aln_sncr_fc_Malig.mdtag.bam -v ~/develop/mutCaller/tests/variants.tsv -t 1 -o test100
~/develop/mutCaller/target/release/mutcaller ALIGNED -b subsampleMalig.sorted.mdtag.bam -v /fh/fast/furlan_s/user/owalt/long_read/lv1/vcf/variants.tsv -t 1 -o out200


~/develop/mutCaller/target/release/mutcaller ALIGNED -b breaksit2.bam -v breaksit.tsv -t 1 -out_long
zcat < ut_long/counts.txt.gz

```

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

First users will generate a simple, 'variants_file' that lists the variants to be counted (tsv.file). The variants file should look something like this:

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
seq     start  ref_nt query_nt name
chr12   112450407   A   G   PTPN11_227A>G
chr12   208248389   G   A   IDH1_132G>A
chr17   7674220 C   T   TP53_248C>T
```

##### Barcodes file

A barcode whitelist.  v3 3prime and v2 5prime 10X whitelists are provided in the data folder.


##### Invocation

**To run mutCaller simply type:**

```plaintext
mutCaller ALIGNED --bam <file.bam> -v <variants.tsv> -o <folder>
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
    -o, --output <output>        output path; will contain mutcaller.log and counts.txt.gz
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
    -o, --output <output>                  output path; will contain mutcaller.log and counts.txt.gz
    -r, --read_len <read_len>              read 2 length (default 90)
    -t, --threads <threads>                threads
    -u, --umi_length <umi_len>             length of umi sequence
    -v, --variants <variants>              path to variants.tsv file (SNVs with mm2 only supported currently) with the
                                           following formatting per line - seqname\tstart\tref_nt\tquery_nt\tname; e.g.
                                           chr12,112450407,A,G,PTPN11_227A>G

```









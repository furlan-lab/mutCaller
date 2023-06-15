# mutCaller
call mutations in 10X data


This is an archived branch of the first version of mutCaller.  This branch is maintained as it retains code for processing 10X data for use with downstream cbsniffer script.  Ultimately it was decided that cbsniffer is too slow, but this code is retained to enable direct head-to-head comparisons of modern mutcaller output with cbsniffer.

To run:

```sh
loc=~/develop/mutCaller # or location where you have cloned the repository
module load SAMtools/1.16.1-GCC-11.2.0 # make samtools accessible on command line
ml R/4.1.2-foss-2021b # make R accessible on command line
cd $loc
git checkout bash
cd mutcaller_rust
git status
cargo build --release
bc=$loc/data/737K-august-2016.txt.gz  #barcode whitelist
out=$loc/cbsniffer
mkdir $out
cd $out
$loc/mutcaller -u -U 10 -b $loc/tests/sequencer_R1.fastq.gz -t $loc/tests/sequencer_R2.fastq.gz -l $bc &&
export fq2=$out/fastq_processed/sample_filtered.fastq &&
export fq3=$out/fastq_processed/sample_filtered_header.fastq &&
$loc/mutcaller_rust/target/release/fastq -t 1 --ifastq ${fq2} > ${fq3}
export transcriptome=/shared/biodata/ngs/Reference/10X/refdata-gex-GRCh38-2020-A #location to cellranger friendly reference
/app/software/CellRanger/6.0.1/lib/bin/STAR --genomeDir $transcriptome/star --readFilesIn ${fq3} --readNameSeparator space --runThreadN 24 --outSAMunmapped Within KeepPairs --outSAMtype BAM SortedByCoordinate &&
#Rscript $loc/scripts/quantReads.R &&
$loc/scripts/addTag.py -u 10 -c 16 Aligned.sortedByCoord.out.bam | samtools view -hbo Aligned.out.tagged.sorted.bam &&
samtools index -@ 24 Aligned.out.tagged.sorted.bam
rm Aligned.sortedByCoord.out.bam &&
rm fastq.log Log.out Log.progress.out SJ.out.tab &&
rm -R fastq_processed &&
rm -R mutcaller &&
rm -R _STARtmp
zcat $bc > $loc/data/737K-august-2016.txt
sed -i 's/$/-1/g' $loc/data/737K-august-2016.txt
cat $loc/data/737K-august-2016.txt | head
python $loc/scripts/cb_sniffer.py Aligned.out.tagged.sorted.bam $loc/tests/variants_cb_sniffer.tsv $loc/data/737K-august-2016.txt test

```

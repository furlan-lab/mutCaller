/**
this module handles the UNALIGNED functions in mutcaller

loc=~/develop/mutCaller
fa=/Users/sfurlan/refs/gencode.v43.transcripts.fa
bc=$loc/data/737K-august-2016.txt.gz 
mutcaller=/Users/sfurlan/develop/mutCaller/target/release/mutcaller
$mutcaller UNALIGNED -v -t 8 -g $fa -b $bc -s $loc/tests/variants.tsv -a bwa -o out_bwa -i $loc/tests/sequencer_R1.fastq.gz -j $loc/tests/sequencer_R2.fastq.gz


bwa mem -t 8 /Users/sfurlan/refs/gencode.v43.transcripts.fa out_bwa/mutcaller_R1.fq.gz > out_bwa/Aligned.sam



#long read JMML woes # make a breaksit file
folder=/fh/scratch/delete90/furlan_s/targ_reseq/230710_longread_JMMLpulldown/1_A01/bc1001--bc1001/1_A01
cd $folder
bam=/fh/scratch/delete90/furlan_s/targ_reseq/230710_longread_JMMLpulldown/1_A01/bc1001--bc1001/1_A01/bc1001--bc1001.mdmapped.bam
variants=/fh/scratch/delete90/furlan_s/targ_reseq/230710_longread_JMMLpulldown/variants.tsv
mutcaller ALIGNED --bam ${bam} -s ${variants} -o mutcaller
samtools view -hb $bam chr1:114705630-114705640 > breaksit_jmml.bam


#copy locally and check after making variants_jmml.tsv
mutcaller ALIGNED --bam ~/develop/mutCaller/tests/breaksit_jmml.bam -s ~/develop/mutCaller/tests/variants_jmml.tsv -o jmml

# grch38
fa=/Users/sfurlan/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa
samtools faidx $fa
bam=~/develop/mutCaller/tests/ds_cellranger.bam
bcftools mpileup -Ou -f $fa $bam | \
bcftools call -Ou -mv > cr.vcf.gz
#fail

cat > cr_var.csv <<EOF
seq,start,ref_nt,query_nt,name
chr1,40861621,A,T,test
EOF
sed 's/,/\t/g' cr_var.csv > cr_var.tsv
cat cr_var.tsv
mutcaller=/Users/sfurlan/develop/mutCaller/target/release/mutcaller
$mutcaller ALIGNED --bam $bam -s cr_var.tsv -o cr -u UB


#try convering to fastq and copying tags if there is a problem with MD tagmaking using calmd
samtools view bc1001--bc1001.corrected.sorted.bam | head
genome=/shared/biodata/ngs/Reference/10x/refdata-gex-GRCh38-2020-A/fasta/genome.fa
samtools fastq -T RG,XC,XM,ac,bc,bq,bx,cx,ec,ls,ma,np,rq,sn,we,ws,zm,qs,qe,XA,rc,gp,nb,CB,CR bc1001--bc1001.corrected.sorted.bam | minimap2 --MD -Y -a $genome -t 10 -o bc1001--bc1001.mapped_mm.sam -
samtools fastq -T RG,XC,XM,ac,bc,bq,bx,cx,ec,ls,ma,np,rq,sn,we,ws,zm,qs,qe,XA,rc,gp,nb,CB,CR breaksit_jmml.bam | minimap2 --MD -Y -a $genome -t 10 -o breaksit_jmml_mapped_mm.sam -
samtools sort breaksit_jmml.sam -O BAM -o breaksit_jmml_sorted.bam
samtools index breaksit_jmml_sorted.bam
rm breaksit_jmml.sam
variants=/fh/scratch/delete90/furlan_s/targ_reseq/230710_longread_JMMLpulldown/variants.tsv
$mutcaller ALIGNED --bam breaksit_jmml_sorted.bam -s ${variants} -o mutcaller
$mutcaller ALIGNED --bam breaksit_jmml_sorted.bam -s variants_jmml.tsv -o mutcaller
## no difference!

**/



#[test]
fn some_test() {
    assert_eq!(2 + 2, 4);
}
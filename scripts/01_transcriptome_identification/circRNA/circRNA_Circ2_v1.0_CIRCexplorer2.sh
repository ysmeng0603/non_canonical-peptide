#!/usr/bin/bash

# $1 accession
# $2 directory stored the unmapped fastq
# $3 directory to store the unmapped sam
# $4 directory to store the circexplore2 parse results and the annotated results
# $5 directory to store the final results，ciri2

index_bwa=./RNA_edit/required_data/bwa_index/hg38
ciri2=./circRNA/ciri_v1/CIRI_v2.0.6/CIRI2.pl
ref_anno=./reference/genome/hg38_ref.txt
fa=./reference/ref_hg38/GRCh38.primary_assembly.genome.fa
gtf=./reference/genome/gencode.hg38.primary_assembly.annotation.gtf

cat $1 | while read id
do
	bwa mem -t 8 -T 19 $index_bwa $2/$id.unmapped.fastq > $3/$id.unmapped.sam 2> $3/$id.bwa.log
	CIRCexplorer2 parse -t BWA -b $4/$id.circexplore2.bed $3/$id.unmapped.sam > $4/$id.parse.log
	CIRCexplorer2 annotate -r $ref_anno -g $fa -b $4/$id.circexplore2.bed -o $4/$id.circexplore2.results.anno
	perl $ciri2 -T 6 -I $3/$id.unmapped.sam -O $5/$id.ciri2_res.txt -F $fa -A $gtf -G $5/$id.ciri2.ciri.log
done

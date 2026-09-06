#!/usr/bin/bash

# $1 sample accession
# $2 whether the fastq is Paired or Single: Paired==P, Single==S
# $3 directory name where stored the fastq.gz files
# $4 dirwctory name where to store the alignment bam files
# $5 directory where to store the anchors.qfa files
# $6 directory where to store the final results
# $7 file contains two columns: Accession       status/disease/cell_lines,etc, tab-separated,example:./circRNA/info/sample_metaInfo.txt


index=./reference/index/bowtie/hg38
chr_fa=./reference/chrom

cat $1 | while read id
do
	if [ $2 == "P" ]
	then
		bowtie2 --sensitive -p 10 -x $index --score-min=C,-15,0 \
			-1 $3/${id}_1.fastq.gz \
                        -2 $3/${id}_2.fastq.gz | samtools view -hbuS -@ 10 - | samtools sort -@ 10 - -o $4/$id.bam

                samtools view -@ 10 -hf 4 $4/$id.bam | samtools view -@ 10 -Sb - >  $4/$id.unmapped.bam
		python ./circbase/unmapped2anchors.py $4/${id}.unmapped.bam > $5/$id.anchors.qfa
                
		bowtie2 --reorder --mm -M20 --score-min=C,-15,0 -x $index -p 10 \
                        -U $5/$id.anchors.qfa | python ./circbase/find_circ.py -r $7 -G $chr_fa -p $id -s $6/$id.sites.log > $6/$id.bed 2> $6/$id.reads

	else

		bowtie2 --sensitive -p 10 -x $index --score-min=C,-15,0 \
			-U $3/$id.fastq.gz | samtools view -@ 10 -hbuS - | samtools sort -@ 10 - -o $4/$id.bam
		 
		samtools view -@ 10 -hf 4 $4/$id.bam | samtools view -@ 10 -Sb - > $4/$id.unmapped.bam
		python ./circbase/unmapped2anchors.py $4/${id}.unmapped.bam > $5/$id.anchors.qfa
		
		bowtie2 --reorder --mm -M20 --score-min=C,-15,0 -x $index -p 10 \
			-U $5/$id.anchors.qfa | python ./circbase/find_circ.py -r $7 -G $chr_fa -p $id -s $6/$id.sites.log > $6/$id.bed 2> $6/$id.reads
	fi
done

#!/usr/bin/bash

# $1 accession
# $2 directory stored the unmapped bam file and to store the out put of star
# $3 the minimum length of the circRNA [int]
# $4 the directory to store the final results

star_index=./reference/index/star
postprocess=./circRNA/GSE89408/circbase/circRNA_finder/postProcessStarAlignment.pl

cat $1 | while read id
do
	bamToFastq -i $2/$id.unmapped.bam -fq $2/$id.unmapped.fastq

	STAR \
		--runThreadN 15 \
		--chimScoreMin 1 \
		--chimSegmentMin 20 \
		--alignIntronMax 100000 \
		--genomeDir $star_index \
		--alignTranscriptsPerReadNmax 100000 \
		--chimOutType Junctions SeparateSAMold \
		--readFilesIn $2/$id.unmapped.fastq \
		--outFilterMismatchNmax 4 \
		--outFilterMultimapNmax 2 \
		--outFileNamePrefix $2/$id
	
	perl $postprocess --starDir $2/$id --minLen $3 --outDir $4
done

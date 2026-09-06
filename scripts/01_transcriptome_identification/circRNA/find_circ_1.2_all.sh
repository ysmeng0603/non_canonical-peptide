#!/usr/bin/bash
#' $1 sample accession
#' $2 directory to store the results
#' $3 directory where to store the anchors.qfa files
#' $4 just as the $7 in the old script

cat $1 | while read id
do
bowtie2 --reorder --mm --score-min=C,-15,0 -x ./circRNA/index/bowti2_index/GRCh38 \
	-q -U $3/$id.anchors.qfa  -p 10 | python ./github/find_circ-master/find_circ.py --genome=./reference/genome/hg38_ref.txt --prefix=$id --reads2samples=$4 --name=Tissue --stats=$2/$id.txt --reads=$2/$id.reads > $2/$id.bed
done

grep CIRCULAR ./venn/venn_need_data1/test/find_circ/$1.bed | \
        grep -v chrM | \
        awk '$5>=2' | \
        grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | \
        ./circRNA/maxlength.py 100000 \
        > ./circRNA/venn/venn_need_data1/test/find_circ/circ_candidates_$1.bed
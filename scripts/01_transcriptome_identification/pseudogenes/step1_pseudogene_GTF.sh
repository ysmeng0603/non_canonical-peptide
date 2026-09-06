#!/bin/bash
cd ./pseudogenes

GTF=./pseudogenes/need_file/gencode.v49.primary_assembly.annotation.gtf
awk '
BEGIN{FS=OFS="\t"}
$0 ~ /^#/ {next}
$9 ~ /gene_type "[^"]*pseudogene"/ {
    print
}
' $GTF > annotation/gencode.v49.pseudogene.clean.gtf
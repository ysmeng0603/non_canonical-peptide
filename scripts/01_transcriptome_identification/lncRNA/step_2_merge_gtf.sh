#!/bin/bash
##
gtf_file="./reference/Homo_sapiens.GRCh38.113.gtf"
merge_list="./lncRNA/StringTie/need_file/mergelist.txt"
#
input_dir="./lncRNA/StringTie/result/step_1_assembly"
out_dir="./lncRNA/StringTie/result/step_2_merge"
mkdir -p $out_dir
#
echo "start:$(date '+%Y-%m-%d %H:%M:%S')"

ls $input_dir/*.gtf > $merge_list

if [[ -s $merge_list ]]; then
    echo "$merge_list exist"
    stringtie --merge \
    -G $gtf_file \
    -o $out_dir/merged.gtf $merge_list
else
    echo "$merge_list does not exist"
fi

echo "finish:$(date '+%Y-%m-%d %H:%M:%S')"

#!/bin/bash

echo RUNNING: generate_fusioninspector_fusions_list.sh cff_file fusion_list_dir
#if [ ! $# -eq 2 ]
#  then
#    echo "incorrect number of arguments supplied"
#    kill $$
#fi
cff_file=$1
fusion_list_file=$2

# $scripts_path/generate_fusioninspector_fusions_list.sh $file $fusion_list_dir
#cff_file_base=$(basename $cff_file)
#awk '{print $14"--"$16 }' $cff_file > $out_dir/$cff_file_base-gene_list.txt 
awk '{print $14"--"$16 }' $cff_file > $fusion_list_file 

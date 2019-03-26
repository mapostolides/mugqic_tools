#!/bin/bash

source setup.sh

scripts_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts
fusion_cluster_file=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/compare_Yue_my_pipelines/RUN_mine/output-Feb26-2019-DEFAULT_PARAMS/fusions/cff/merged.cff.reann.dnasupp.bwafilter.30.cluster

output_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts/TEST/generate_category_table/output
mkdir -p $output_dir

python $scripts_dir/generate_category_table.py $fusion_cluster_file $output_dir 

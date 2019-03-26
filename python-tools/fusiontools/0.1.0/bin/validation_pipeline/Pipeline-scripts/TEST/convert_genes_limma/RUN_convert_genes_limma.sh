#!/bin/bash

scripts_path=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts
cff_file=SJHGG136_D.integrate.cff
test_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts/TEST/convert_genes_limma
$scripts_path/rename_cff_file_genes.py $cff_file > $test_dir/SJHGG136_D.integrate.cff.renamed

#/hpf/tools/centos6/R/3.5.1/bin/Rscript $code_dir/convert_genes_limma.R $cff_file

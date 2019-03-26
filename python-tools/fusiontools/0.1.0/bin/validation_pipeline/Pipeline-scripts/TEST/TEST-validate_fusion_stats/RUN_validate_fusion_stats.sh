#!/bin/bash

source setup.sh 

scripts_path=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts

validation_pipeline_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline

validated_fusions_file=$validation_pipeline_dir/sim45_validation/sim45_validated_fusions_infofile-GENE_NAMES-bedfile_intersected.tsv 

cluster_file_bed_annotation=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/smc_rna_sim45/output-Feb15-2019-NO_fusioninspector_ericscript_param_mod/fusions/cff/merged.cff.reann.dnasupp.bwafilter.30.cluster.bed_gene_names-TWEAK1

pipeline_outdir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts/TEST/TEST-validate_fusion_stats/outdir

python $scripts_path/validate_fusion_stats.py $validated_fusions_file $cluster_file_bed_annotation $pipeline_outdir 

#!/bin/bash
source setup.sh
scripts_path=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts
#python $scripts_path/annotate_called_fusion_file.py $pipeline_outdir/$basename_cff.renamed.cluster $sample_name $pipeline_outdir/$basename_cff.renamed.bed_gene_names.cluster

#TEST ANNOTATION OF CFF FILE
python $scripts_path/CFF_INPUT-annotate_called_fusion_file.py inversion.cff test_sample ./outfiles/inversion.bed_gene_names.cff

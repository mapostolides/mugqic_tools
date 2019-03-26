#!/bin/bash
source setup.sh
scripts_path=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts
#python $scripts_path/merge_called_fusions.py $pipeline_outdir/$basename_cff.renamed.bed_gene_names.cluster $sample_name $pipeline_outdir/$basename_cff.renamed.bed_gene_names.merged_breakpoints.cluster

#TEST TO SEE THAT INVERTED GENE PAIRS GET MERGED (i.e. GeneA-GeneB and GeneB-GeneA are merged)
python $scripts_path/merge_called_fusions.py st-jude-STAR_FUSION.cff.renamed.bed_gene_names.cluster st-jude st-jude-STAR_FUSION.cff.renamed.bed_gene_names.merged_breakpoints.cluster

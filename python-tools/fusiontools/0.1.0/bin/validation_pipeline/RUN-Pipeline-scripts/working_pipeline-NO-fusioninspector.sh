#!/bin/bash

#bash test_script.sh $sample_name $cff_file $validation_pipeline_dir $pipeline_outdir $validated_fusions_file 
#CONFIGURATION
# run setup in "scripts" directory of my home folder
source /hpf/largeprojects/ccmbio/mapostolides/scripts/setup.sh
sample_name=$1
cff_file=$2
validation_pipeline_dir=$3
pipeline_outdir=$4
validated_fusions_file=$5
scripts_path=$validation_pipeline_dir/Pipeline-scripts

#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts/rename_cff_file_genes.py /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/fusion_pipeline_output_files/$basename_cff

#$STAR_Fusion_outfile_path=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/STAR-Fusion/STAR-Fusion_outputDec13-v1.5.0-skewer_trimmed_reads/star-fusion.fusion_predictions.abridged.tsv

# rename genes, write new cff file
# original input file: testfiles/$sample_name-tune_run_1_merged.cff
basename_cff=$(basename $cff_file)
#$scripts_path/rename_cff_file_genes.py $cff_file > $pipeline_outdir/$basename_cff.renamed

# swap columns to create valid input file for validate_fusion_stats.py (CategoryFusions line format)
#$scripts_path/convert_cff_to_fake_cluster.sh $pipeline_outdir/$basename_cff.renamed > $pipeline_outdir/$basename_cff.renamed.cluster

# intersect called fusion breakpoints with .bed file located here:  "/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID.bed"
# add 4 columns to .cluster file: gene1_candidates        gene2_candidates        gene1_strands   gene2_strands,
# where "candidates" refers to the possibility that a breakpoint will intersect multiple genes on either strand, and it may not be clear which gene is participating in the fusion transcript
#if [ ! -f $pipeline_outdir/$basename_cff.renamed.bed_gene_names.cluster ];
#    then
#    python $scripts_path/CLUSTER_INPUT-annotate_called_fusion_file.py $pipeline_outdir/$basename_cff.renamed.cluster $sample_name $pipeline_outdir/$basename_cff.renamed.bed_gene_names.cluster 
#else
#    echo "   ANNOTATED BREAKPOINT FILE EXISTS"
#fi

python $scripts_path/merge_called_fusions.py $pipeline_outdir/$basename_cff.renamed.bed_gene_names.cluster $sample_name $pipeline_outdir/$basename_cff.renamed.bed_gene_names.merged_breakpoints.cluster   


# run validation script on .cluster file
python $scripts_path/validate_fusion_stats.py $validated_fusions_file $pipeline_outdir/$basename_cff.renamed.bed_gene_names.merged_breakpoints.cluster $sample_name $pipeline_outdir



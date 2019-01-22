#!/bin/bash

#CONFIGURATION
# run setup in "scripts" directory of my home folder
source setup.sh
sample_name=smc_rna_sim45

#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts/rename_cff_file_genes.py /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/fusion_pipeline_output_files/$sample_name-merged_star_fusion.cff
DATE=`date +%Y-%m-%d`
validation_pipeline_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline
pipeline_outdir=$validation_pipeline_dir/validation_pipeline_output_dirs/$sample_name-NO_FusionInspector-$DATE
mkdir $pipeline_outdir 


scripts_path=$validation_pipeline_dir/Pipeline-scripts

STAR_Fusion_outfile_path=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/STAR-Fusion/STAR-Fusion_outputDec13-v1.5.0-skewer_trimmed_reads/star-fusion.fusion_predictions.abridged.tsv
python $scripts_path/convert_STAR-Fusion_output_to_cff.py $STAR_Fusion_outfile_path $pipeline_outdir/$sample_name-star_fusion.cff

# merge fusion pipeline output cff with star_fusion cff
#cat $validation_pipeline_dir/fusion_pipeline_output_files/$sample_name-merged.cff $pipeline_outdir/$sample_name-star_fusion.cff > $pipeline_outdir/$sample_name-merged_star_fusion.cff
#cat $validation_pipeline_dir/fusion_pipeline_output_files/$sample_name-merged-no-ericscript-integrate.cff $pipeline_outdir/star_fusion.cff > $pipeline_outdir/$sample_name-merged_star_fusion.cff

# rename genes, write new cff file
# original input file: testfiles/$sample_name-tune_run_1_merged.cff
$scripts_path/rename_cff_file_genes.py $pipeline_outdir/$sample_name-merged_star_fusion.cff > $pipeline_outdir/$sample_name-merged_star_fusion.cff.renamed

# swap columns to create valid input file for validate_fusion_stats.py (CategoryFusions line format)
$scripts_path/convert_cff_to_fake_cluster.sh $pipeline_outdir/$sample_name-merged_star_fusion.cff.renamed > $pipeline_outdir/$sample_name-merged_star_fusion.cff.renamed.cluster

# intersect called fusion breakpoints with .bed file located here:  "/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID.bed"
# add 4 columns to .cluster file: gene1_candidates        gene2_candidates        gene1_strands   gene2_strands,
# where "candidates" refers to the possibility that a breakpoint will intersect multiple genes on either strand, and it may not be clear which gene is participating in the fusion transcript
if [ ! -f $pipeline_outdir/$sample_name-merged_star_fusion.cff.renamed.bed_gene_names.cluster ];
    then
    python $scripts_path/annotate_called_fusion_file.py $pipeline_outdir/$sample_name-merged_star_fusion.cff.renamed.cluster $sample_name $pipeline_outdir/$sample_name-merged_star_fusion.cff.renamed.bed_gene_names.cluster 
else
    echo "   ANNOTATED BREAKPOINT FILE EXISTS"
fi

python $scripts_path/merge_called_fusions_by_breakpoint.py $pipeline_outdir/$sample_name-merged_star_fusion.cff.renamed.bed_gene_names.cluster $sample_name $pipeline_outdir/$sample_name-merged_star_fusion.cff.renamed.bed_gene_names.merged_breakpoints.cluster   

# run validation script on .cluster file
python $scripts_path/validate_fusion_stats.py $validation_pipeline_dir/sim45_validated_fusions_infofile-GENE_NAMES.tsv $pipeline_outdir/$sample_name-merged_star_fusion.cff.renamed.bed_gene_names.merged_breakpoints.cluster $sample_name $pipeline_outdir



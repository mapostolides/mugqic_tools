#!/bin/bash

#$pipeline_rundir/working_pipeline-NO-fusioninspector.sh $sample_name $cff_file $validation_pipeline_dir $pipeline_outdir $validated_fusions_file 
sample_name=$1
cff_file=$2
validation_pipeline_dir=$3
pipeline_outdir=$4
validated_fusions_file=$5
scripts_path=$validation_pipeline_dir/Pipeline-scripts

#SETUP
source setup.sh # run setup in "scripts" directory of my home folder

renamed_dir=$pipeline_outdir/$sample_name-callers-separate-renamed
mkdir -p $renamed_dir 

fusion_list_dir=$pipeline_outdir/fusioninspector_fusion_list_files
mkdir -p $fusion_list_dir

FI_filtered_dir=$pipeline_outdir/FI-filtered_cff_files
mkdir -p $FI_filtered_dir

# CONFIGURATION
pre_fusioninspector=0
run_fusioninspector=0
filter_fusioninspector_output=0
cat_FI_cff_files=1


left_fastq=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/smc_rna_sim45/smc_rna_sim45.pair1.fastq.gz                                                                       
right_fastq=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/smc_rna_sim45/smc_rna_sim45.pair2.fastq.gz  


#STAR_Fusion_outfile_path="/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/STAR-Fusion/STAR-Fusion_outputDec6-v1.5.0/star-fusion.fusion_predictions.abridged.tsv"
#STAR_Fusion_outfile_path="/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/STAR-Fusion/STAR-Fusion_outputDec13-v1.5.0-skewer_trimmed_reads/star-fusion.fusion_predictions.abridged.tsv"

#BEGIN pre_fusioninspector
if [ $pre_fusioninspector -eq 1 ]
then
    echo running pre_fusioninspector steps

    #GENERATE STAR-FUSION CFF
    #echo generate STAR-Fusion .cff file
    #star_fusion_outfile=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/STAR-Fusion/STAR-Fusion_outputDec13-v1.5.0-skewer_trimmed_reads/star-fusion.fusion_predictions.abridged.tsv
    #star_fusion_cff_file=$validation_pipeline_dir/fusion_pipeline_output_files/$sample_name-callers-separate/$sample_name-star_fusion.cff
    #python $scripts_path/convert_STAR-Fusion_output_to_cff.py $star_fusion_outfile $star_fusion_cff_file
    
    
    #RENAME CALLED FUSION CFF FILES FOR ALL CALLERS SEPARATELY
    echo generate renamed fusion gene .cff files for all callers
    #separate_callers_dir already exists in fusion_pipeline_output_files directory
    separate_callers_dir=$validation_pipeline_dir/fusion_pipeline_output_files/$sample_name-callers-separate
    for file in $separate_callers_dir/*;do
         file=$(basename $file)
         file_renamed="${file%%.*}".renamed.cff
         $scripts_path/rename_cff_file_genes.py $separate_callers_dir/$file > $renamed_dir/$file_renamed 
    done
    
    #GENERATE FUSIONINSPECTOR FUSION GENE PAIR LIST FILES
    echo generate FusionInspector fusion gene pair list files for all callers
    for cff_file in $renamed_dir/*.renamed.cff;do
        $scripts_path/generate_fusioninspector_fusions_list.sh $cff_file $fusion_list_dir 
    done

#END pre_fusioninspector
fi

# SUBMIT FUSIONINSPECTOR JOB FOR EACH CALLER
FI_outdir=$pipeline_outdir/$sample_name-fusioninspector_output
mkdir -p $FI_outdir
FI_errdir=$pipeline_outdir/$sample_name-fusioninspector_ERROR
mkdir -p $FI_errdir

if [ $run_fusioninspector -eq 1 ]
then
    echo $fusion_list_dir
    #ls $fusion_list_dir 
    for fusion_list_file in $fusion_list_dir/*; do 
        #echo $fusion_list_file
        #fusion_list_basename=$(basename $fusion_list_file)
        #echo $fusion_list_basename
        #caller_name="${fusion_list_basename%%.*}"
        #echo $caller_name
        JOB=$(qsub -v left_fastq=$left_fastq,right_fastq=$right_fastq,fusion_list_file=$fusion_list_file,pipeline_outdir=$pipeline_outdir,FI_outdir=$FI_outdir -o $FI_errdir -e $FI_errdir $scripts_path/run_fusioninspector-generic-qsub.sh)
        echo $JOB 
    done
fi

# FILER FUSONS USING FUSIONINSPECTOR OUTPUT


if [ $filter_fusioninspector_output -eq 1 ]
    then

    # DEFINE SUMMARY FILE PATH
    FI_summary_file=$FI_filtered_dir/$sample_name-FI-summary.txt
    rm $FI_summary_file
    printf caller_name'\t'all_caller_calls'\t'FI_total_calls'\t'FI_validated_cff_fusions'\t'unvalidated_cff_fusions'\n' >> $FI_summary_file

    #FILTER
    for caller_name in defuse ericscript fusionmap integrate star_fusion; do
        caller_outdir=$FI_outdir/$sample_name-$caller_name-FI_output
        FI_output_file=$caller_outdir/finspector.fusion_predictions.final.abridged.FFPM
        renamed_cff_file=$renamed_dir/$sample_name-$caller_name\.renamed.cff
        
        # DEFINE UNVALIDATED AND VALIDATED CFF FILE PATHS
        file=$(basename $renamed_cff_file)
        FI_validated="${file%%.*}".renamed.FI-validated.cff
        FI_validated_cff_file=$FI_filtered_dir/$FI_validated
        FI_unvalidated="${file%%.*}".renamed.FI-unvalidated.cff
        FI_unvalidated_cff_file=$FI_filtered_dir/$FI_unvalidated
        
        python $scripts_path/filter_fusions_using_fusioninspector_output.py $FI_output_file $renamed_cff_file $FI_validated_cff_file $FI_unvalidated_cff_file $FI_summary_file $caller_name
    done

    # ADD DESCRIPTIONS TO FI SUMMARY FILE
    echo >> $FI_summary_file 
    echo DESCRIPTIONS OF COLUMN HEADERS >> $FI_summary_file 
    echo all_caller_calls: Total number of calls made by fusion caller, before filtering by FusionInspector >> $FI_summary_file
    echo FI_total_calls: The number of called fusions that FusionInspector detected. This is based on input gene pair list>> $FI_summary_file
    echo FI_validated_cff_fusions: The number of fusions called by a given fusion caller that are categorized as valid by FusionInspector >> $FI_summary_file  
    echo unvalidated_cff_fusions: The number of called fusions which are NOT detected as valid fusions by FusionInspector>> $FI_summary_file  
fi

# MERGE ALL FUSIONINSPECTOR FILTERED FILES INTO ONE CFF FILE
if [ $cat_FI_cff_files -eq 1 ]
    then
    defuse=$FI_filtered_dir/$sample_name-defuse.renamed.FI-validated.cff
    ericscript=$FI_filtered_dir/$sample_name-ericscript.renamed.FI-validated.cff
    fusionmap=$FI_filtered_dir/$sample_name-fusionmap.renamed.FI-validated.cff
    integrate=$FI_filtered_dir/$sample_name-integrate.renamed.FI-validated.cff
    star_fusion=$FI_filtered_dir/$sample_name-star_fusion.renamed.FI-validated.cff
    cat $defuse $ericscript $fusionmap $integrate $star_fusion > $FI_filtered_dir/total.renamed.FI-validated.cff
#    cat defuse ericscript fusionmap integrate star_fusion; do
fi

# swap columns to create valid input file for validate_fusion_stats.py (CategoryFusions line format)
$scripts_path/convert_cff_to_fake_cluster.sh $FI_filtered_dir/total.renamed.FI-validated.cff > $FI_filtered_dir/total.renamed.FI-validated.cff.cluster

# intersect called fusion breakpoints with .bed file located here:  "/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID.bed"
# add 4 columns to .cluster file: gene1_candidates        gene2_candidates        gene1_strands   gene2_strands,
# where "candidates" refers to the possibility that a breakpoint will intersect multiple genes on either strand, and it may not be clear which gene is participating in the fusion transcript
if [ ! -f $pipeline_outdir/$sample_name\.total.cff.renamed.FI-validated.bed_gene_names.cluster ];
    then
    python $scripts_path/annotate_called_fusion_file.py $FI_filtered_dir/total.renamed.FI-validated.cff.cluster smc_rna_sim45 $pipeline_outdir/$sample_name\.total.cff.renamed.FI-validated.bed_gene_names.cluster 
else
    echo "   ANNOTATED FILE EXISTS"
fi

python $scripts_path/merge_called_fusions_by_breakpoint.py $pipeline_outdir/$sample_name\.total.cff.renamed.FI-validated.bed_gene_names.cluster smc_rna_sim45 $pipeline_outdir/$sample_name\.total.cff.renamed.FI-validated.bed_gene_names.merged_breakpoints.cluster 

# run validation script on .cluster file
python $scripts_path/validate_fusion_stats.py $validation_pipeline_dir/sim45_validated_fusions_infofile-GENE_NAMES.tsv $pipeline_outdir/$sample_name\.total.cff.renamed.FI-validated.bed_gene_names.merged_breakpoints.cluster  smc_rna_sim45 $pipeline_outdir


# SCRAP
#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts/rename_cff_file_genes.py /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/fusion_pipeline_output_files/smc_rna_sim45-merged_star_fusion.cff

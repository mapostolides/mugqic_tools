#!/bin/bash
#PBS -N FusionInspector
#PBS -l nodes=1:ppn=4
#PBS -l mem=60g
#PBS -l vmem=60g
#PBS -l walltime=24:00:00


module load samtools
#module load trinityrnaseq
module load python/2.7.11
module load gcc/5.2.0
module load perl/5.28.0_mikeapos

export PATH=$PATH:/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/FusionInspector-run/FusionInspector # FusionInspector
export PATH=$PATH:/hpf/largeprojects/ccmbio/mapostolides/MODULES/bin #gmap-2018-07-04, htslib-1.3 (for GMAP and bgzip dependencies, respectively)
export PATH=$PATH:/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/FusionInspector-run/dependencies/trinityrnaseq-Trinity-v2.8.4
export PATH="/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/STAR/STAR-2.6.1c/source:$PATH"

ctat_lib_dir=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/GRCh37_v19_CTAT_lib_Feb092018.plug-n-play/ctat_genome_lib_build_dir
#top_dir=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/FusionInspector-run

# variables from qsub -v: pipeline_outdir,FI_outdir,left_fastq,right_fastq,fusion_list_file

# CREATE SEPARATE OUTPUT DIRECTORY FOR CALLER
fusion_list_basename=$(basename $fusion_list_file)
caller_name="${fusion_list_basename%%.*}" 
caller_outdir=$FI_outdir/$caller_name-FI_output
mkdir -p $caller_outdir

FusionInspector --fusions $fusion_list_file \
                --genome_lib $ctat_lib_dir \
                --left_fq $left_fastq --right_fq $right_fastq \
                --out_dir $caller_outdir \
                --out_prefix finspector \
                --prep_for_IGV

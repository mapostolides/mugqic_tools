#!/bin/bash
#NOTES
#$reann_test_dir/cff_files/NOTES-merged.cff_RPRD2--LAMC2

#CONSTANTS
#gene_info_file=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/Homo_sapiens.gene_info
info_file_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/rename_genes_test/info_file_subsets
reann_test_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing
test_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline
outdir=$reann_test_dir/outdir-April-30-2020
#gene_bed_total=$reann_test_dir/ens_known_genes.renamed.bed
gene_bed_total=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/annotation/ens_known_genes.renamed.ENSG.bed
bed_dir=$reann_test_dir/bed_files
cff_dir=$reann_test_dir/cff_files
#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion.py 
#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/generate_common_fusion_stats.py 
run_pipeline (){
  outdir=$1
  mkdir -p $outdir
  cff=$2
  gene_bed=$3
  gene_info_file=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/Homo_sapiens.gene_info
  truth_fusions=$4
  fusiontools=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin
  genome_fasta=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/human_g1k_v37_decoy.fasta

  #Rename cff
  #python $fusiontools/rename_cff_file_genes-GENAP.py $cff $gene_info_file > $outdir/$(basename $cff).renamed
  cff=$outdir/$(basename $cff).renamed

  #Annotate cff
  #python $fusiontools/reann_cff_fusion.py $cff $gene_bed $genome_fasta > $outdir/$(basename $cff).reann 
  cff=$outdir/$(basename $cff).reann

  #Merge
  cluster=$outdir/$(basename $cff).cluster
  #sh $fusiontools/RUN_cluster_genes_breakpoints.sh $cff $outdir > $cluster

  #Benchmark
  #/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster.sh $outdir $truth_fusions $cff $cluster true true true true false
  /hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster.sh $outdir $truth_fusions $cff $cluster false true true true false
}

date=June-1-2020

# DIPG -- benchmark experimentally validated in RNA
outdir=$test_dir/DIPG.RNA_validated.benchmark.$date
gene_bed=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/annotation/ens_known_genes.renamed.ENSG.bed
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/DIPG.RNA_validated.renamed.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/DIPG/run_DIPG_samples/DIPG_output_7CALLERS_April-17-2020/fusions/cff/merged.cff
#run_pipeline $outdir $cff $gene_bed $truth_fusions

#BT474.KPL4.MCF7.SKBR3
outdir=$test_dir/BT474.KPL4.MCF7.SKBR3.benchmark.$date
gene_bed=$gene_bed_total
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/BT474.KPL4.MCF7.SKBR3.truth_set.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output.BT474.KPL4.MCF7.SKBR3-April-9-2020/fusions/cff/merged.cff
run_pipeline $outdir $cff $gene_bed $truth_fusions
exit 0

# SIM50 2500 fusions files:
outdir=$test_dir/SIM50.2500_TP.benchmark.$date
gene_bed=$gene_bed_total
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/sim50.truth_set.2500_fusions.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM-April-17-2020/fusions/cff/merged.cff
run_pipeline $outdir $cff $gene_bed $truth_fusions

#SIM101 2500 fusions files, same truth set as SIM50
outdir=$test_dir/SIM101.2500_TP.benchmark.$date
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/sim_101.truth_set.2500_fusions.dat
gene_bed=$gene_bed_total
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM_101-April-23-2020/fusions/cff/merged.cff
run_pipeline $outdir $cff $gene_bed $truth_fusions

#NEGATIVE CONTROL BEERS
outdir=$test_dir/NEG_CONTROL_BEERS.benchmark.$date
gene_bed=$gene_bed_total
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/negative_control_beers/output-7CALLERS-April-17-2020/fusions/cff/merged.cff
run_pipeline $outdir $cff $gene_bed $truth_fusions

# SIM45.SIM52.combined
gene_bed=$gene_bed_total
outdir=$test_dir/SIM45.SIM52.benchmark.$date
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/sim45.sim52.combined/output-7CALLERS-April-22-2020/fusions/cff/merged.cff
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/renamed_truth_sets/sim45.sim52.truth_set.dat
run_pipeline $outdir $cff $gene_bed $truth_fusions


exit 0
#SIM50 benchmark only- test using null truth set
outdir=$test_dir/SIM50.2500.benchmark_null_truth.May-26-2020
truth_fusions=$outdir/NULL_TRUTH.dat
cff_orig=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/SIM50.2500.benchmark.May-25-2020_genebed_renamed.ENSG_gene_bed/merged.cff.renamed.reann
cluster=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/SIM50.2500.benchmark.May-25-2020_genebed_renamed.ENSG_gene_bed/merged.cff.renamed.reann.cluster
/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster.sh $outdir $truth_fusions $cff_orig $cluster true true true true false
exit 0


# DIPG
outdir=$test_dir/DIPG.May-26-2020_genebed_renamed.ENSG_gene_bed
gene_bed=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/annotation/ens_known_genes.renamed.ENSG.bed
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/haas_2019_simulated_dataset_FILES/sim_50.truth_set.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/DIPG/run_DIPG_samples/DIPG_output_7CALLERS_April-17-2020/fusions/cff/merged.cff
#run_pipeline $outdir $cff $gene_bed $truth_fusions
#exit 0

# SIM50 2500 fusions files:
#ENSG gene_bed
outdir=$test_dir/SIM50.2500.benchmark.May-25-2020_genebed_renamed.ENSG_gene_bed.pygeneann_April_27_2020
gene_bed=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/annotation/ens_known_genes.renamed.ENSG.bed
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/haas_2019_simulated_dataset_FILES/sim_50.truth_set.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM-April-17-2020/fusions/cff/merged.cff
run_pipeline $outdir $cff $gene_bed $truth_fusions
exit 0
#ALL gene_bed
outdir=$test_dir/SIM50.2500.benchmark.May-25-2020_genebed_renamed.ALL_gene_bed
gene_bed=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/annotation/ens_known_genes.renamed.bed
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/haas_2019_simulated_dataset_FILES/sim_50.truth_set.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM-April-17-2020/fusions/cff/merged.cff
run_pipeline $outdir $cff $gene_bed $truth_fusions


#SIM45.SIM52
gene_bed=$gene_bed_total
outdir=$test_dir/SIM45.SIM52.benchmark.May-25-2020_genebed_renamed.ENSG
truth_fusions=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/sim45.sim52.combined/annotate_TP_set/sim45.sim52.truth_set.KEEP.dat
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/sim45.sim52.combined/output-7CALLERS-April-22-2020/fusions/cff/merged.cff
run_pipeline $outdir $cff $gene_bed $truth_fusions
exit 0


outdir=$test_dir/DIPG.TEST_FUSIONS.May-20-2020
mkdir $outdir
gene_bed=$gene_bed_total
cff=/hpf/largeprojects/ccmbio/mapostolides/DIPG/run_DIPG_samples/DIPG_output_7CALLERS_April-17-2020/fusions/cff/merged.cff
#run_pipeline $outdir $cff $gene_bed $truth_fusions

outdir=$test_dir/outdir.TEST_FUSIONS.May-20-2020
#BCL3--CTB-171A8.1
#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files/BCL3--CTB-171A8.1.cff
cff=$cff_dir/BCL3--CTB-171A8.1.cff
gene_bed=$bed_dir/BCL3--CTB-171A8.1.bed
info_file=$info_file_dir/Homo_sapiens.gene_info.BCL3
#run_pipeline $outdir $cff $gene_bed


echo ARRIBA
cff=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/rename_genes_test/arriba.weird_brackets.cff
#python $fusiontools/rename_cff_file_genes-GENAP.py $cff $gene_info_file
run_pipeline $outdir $cff $gene_bed
echo INTEGRATE
cff=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/rename_genes_test/integrate_slash.delimeter.cff
run_pipeline $outdir $cff $gene_bed
#python $fusiontools/rename_cff_file_genes-GENAP.py $cff $gene_info_file

# TEST FUSIONS
#HES1--IL29
cff=$reann_test_dir/HES1--IL29/HES1--IL29.cff
gene_bed=$reann_test_dir/ens_known_genes_HES1_IL29.bed
run_pipeline $outdir $cff $gene_bed
#SAA1--ATP1A1
cff=$reann_test_dir/SAA1--ATP1A1/SAA1--ATP1A1.cff
gene_bed=$reann_test_dir/ens_known_genes_SAA1--ATP1A1.bed
run_pipeline $outdir $cff $gene_bed
#BCL3--CTB-171A8.1
cff=$cff_dir/BCL3--CTB-171A8.1.cff
gene_bed=$bed_dir/BCL3--CTB-171A8.1.bed
run_pipeline $outdir $cff $gene_bed
#ANKIB1--AKAP9
cff=$reann_test_dir/cff_files/merged.cff_AKAP9--ANKIB1
gene_bed=$reann_test_dir/ens_known_genes.AKAP9--ANKIB1.bed
run_pipeline $outdir $cff $gene_bed
#FOSB--AADACL2
cff=$reann_test_dir/FOSB--AADACL2/FOSB--AADACL2.cff
gene_bed=$bed_dir/ens_known_genes.MIR548H2.FOSB.AADACL2.bed
run_pipeline $outdir $cff $gene_bed

exit 0

#SAA1--ATP1A1
outdir=$reann_test_dir/SAA1--ATP1A1
cff=$reann_test_dir/SAA1--ATP1A1/SAA1--ATP1A1.cff
#run_pipeline $outdir $cff $gene_bed_total

#SIM45.SIM52
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/sim45.sim52.combined/output-7CALLERS-April-22-2020/fusions/cff/merged.cff
#outdir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/TEST_bwafilter-May-4-2020/SIM45.SIM52.cff_10000_bp.May-6-2020
outdir=$reann_test_dir/SIM45.SIM52.May-12-2020
run_pipeline $outdir $cff $gene_bed_total

#BEERS -ive control
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/montreal_benchmark_validation/negative_control_beers/output-7CALLERS-April-17-2020/fusions/cff/merged.cff
outdir=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/neg.control.beers.May-5-2020
#run_pipeline $outdir $cff $gene_bed_total

##BT474.KPL4.MCF7.SKBR3
outdir=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/BT474.KPL4.MCF7.SKBR3.cff-April-30-2020
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output.BT474.KPL4.MCF7.SKBR3-April-9-2020/fusions/cff/merged.cff
#run_pipeline $outdir $cff $gene_bed_total

##SIM50
outdir=/hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/SIM50.cff-April-30-2020
cff=/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/haas_2019_validation/output-SIM-April-17-2020/fusions/cff/merged.cff
#run_pipeline $outdir $cff $gene_bed_total

##
cff=$reann_test_dir/cff_files/merged.cff_BCL3 

#BCL3--CTB-171A8.1
cff=$cff_dir/BCL3--CTB-171A8.1.cff
gene_bed=$bed_dir/BCL3--CTB-171A8.1.bed
gene_bed=$gene_bed_total
#ANKIB1--AKAP9
#gene_bed=$reann_test_dir/ens_known_genes.ANKIB1--AKAP9.bed
#cff=$reann_test_dir/cff_files/merged.cff_AKAP9--ANKIB1
#run_pipeline $outdir $cff $gene_bed

#FOSB--AADACL2
cff=$reann_test_dir/FOSB--AADACL2/FOSB--AADACL2.cff
gene_bed=$bed_dir/ens_known_genes.MIR548H2.FOSB.AADACL2.bed
#1 entry, failing due to breakpt being 1 bp off:
#cff=$reann_test_dir/FOSB--AADACL2/FOSB--AADACL2.cff.test
#run_pipeline $outdir $cff $gene_bed


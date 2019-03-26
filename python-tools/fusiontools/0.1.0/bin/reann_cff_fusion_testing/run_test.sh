#!/usr/bin/env bash
echo USAGE: ../run_test.sh \$1=\*.cff \$2=ens_known_genes.bed human_g1k_v37_decoy.fasta
#python -m trace --listfuncs ../reann_cff_fusion.py $1 $2 human_g1k_v37_decoy.fasta > test 
#python -m trace --listfuncs ../reann_cff_fusion.py $1 $2 human_g1k_v37_decoy.fasta

known_genes_file_name=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID.bed

/Users/mapostolides/miniconda3/envs/pysam/bin/python ../reann_cff_fusion.py $1  human_g1k_v37_decoy.fasta

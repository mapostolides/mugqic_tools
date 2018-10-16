source setup.sh
echo USAGE: ../run_test.sh \$1=\*.cff \$2=ens_known_genes.bed human_g1k_v37_decoy.fasta
#python -m trace --listfuncs ../reann_cff_fusion.py $1 $2 human_g1k_v37_decoy.fasta > test 
python -m trace --listfuncs ../reann_cff_fusion.py $1 $2 human_g1k_v37_decoy.fasta

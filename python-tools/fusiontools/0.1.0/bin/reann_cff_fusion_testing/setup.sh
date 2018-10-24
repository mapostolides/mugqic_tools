#!/bin/bash

# Set up environment for running MUGQIC pipelines:
# 2.2.0 is neded for fusion pipeline because it imports "pysam". latest_hpf does not work for this reason
module load mugqic-pipelines/2.2.0
#module load mugqic-pipelines
#export PYTHONPATH=/hpf/largeprojects/ccmbio/smark/mugqic_tutorial/:$PYTHONPATH
export PYTHONPATH=/hpf/largeprojects/ccmbio/mapostolides/LIS_fusion/:$PYTHONPATH
export MODULEPATH=$MODULEPATH:/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/modulefiles
module load /hpf/largeprojects/ccmbio/mapostolides/gene_fusion/modulefiles/fusiontools/0.1.0 
#module load fusiontools

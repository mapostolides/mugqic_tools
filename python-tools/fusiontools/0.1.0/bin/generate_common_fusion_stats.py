#!/usr/bin/env python
import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/Genap_ccm/pygenefusionann/")
#import pygeneann
import pygeneann_OLD_star_fusion_defuse_style as pygeneann
#import pygeneann_STEPH as pygeneann
import sequtils
import argparse


testing=0
if testing:
    cff_file = "/Users/mapostolides/Desktop/mugqic_tools/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/cff_files_May6/outdir/ANKIB1--AKAP9.cff.reann"
else:
    parser = argparse.ArgumentParser()

    parser.add_argument('cff_file', action='store', help='CFF file')

    args = parser.parse_args()
    cff_file=args.cff_file

cffstats = pygeneann.CffFusionStats(cff_file)

cffstats.generate_common_fusion_stats_by_genes(cff_file)


#!/usr/bin/env python
import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/Genap_ccm/pygenefusionann/")
#import pygeneann
import pygeneann_OLD_star_fusion_defuse_style as pygeneann
import sequtils
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('cff_file', action='store', help='CFF file')

args = parser.parse_args()

cffstats = pygeneann.CffFusionStats(args.cff_file)	

cffstats.generate_common_fusion_stats_by_genes(args.cff_file)


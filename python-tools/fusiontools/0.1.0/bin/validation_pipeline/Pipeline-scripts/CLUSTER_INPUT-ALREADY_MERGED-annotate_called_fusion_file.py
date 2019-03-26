#! /usr/bin/env python

"""
Add the following to called fusions file: gene names, category of gene, and region of gene within which breakpoints are located
"""

import sys
import argparse
from validation_classes import *
from pygeneann import *
parser = argparse.ArgumentParser()

parser.add_argument('cluster_file', action='store', help='A file containing a set of called fusions')
parser.add_argument('annotation_file', action='store', help='known_genes_file_name')
parser.add_argument('output_file', action='store', help='absolute path of the output file')  

args = parser.parse_args()


def map_called_fusions_to_known_genes(called_fusions_file_name, gene_ann_bed):
    """
    Maps start position of called fusions to a known gene in the known genes .bed file
    """
    genebed_objects = []
    called_fusions = []
    for line in open(gene_ann_bed, 'r+'):
        genebed_object = GeneBed(line)
        genebed_objects.append(genebed_object)
    for line in open(called_fusions_file_name):
        if line.startswith('#'): continue
        fusion = CategoryFusions(line)
        called_fusions.append(fusion)

    print("loaded called fusions and GeneBed objects")
    gene_pairs = []
    for fusion in called_fusions:
        # check whole .bed file for candidade gene mappings
        for genebed in genebed_objects:
            if ( genebed.start <= fusion.breakpoint_1[0] <= genebed.end  ) and (fusion.chr1 == genebed.chr):
                #print(genebed.start,  fusion.breakpoint_1[0], genebed.end )
                fusion_left_names = [tup[0] for tup in fusion.left]
                if genebed.gene_name not in fusion_left_names:
                    fusion.left.append((genebed.gene_name, genebed.type, genebed.strand, [genebed.start, genebed.end]))
            if (  genebed.start <= fusion.breakpoint_2[0] <= genebed.end)  and (fusion.chr2 == genebed.chr):
                fusion_right_names = [tup[0] for tup in fusion.right]
                if genebed.gene_name not in fusion_right_names:
                    fusion.right.append((genebed.gene_name, genebed.type, genebed.strand, [genebed.start, genebed.end]))
        gene_pairs.append((fusion.left, fusion.right))
        print(fusion.left, fusion.right)
    return (gene_pairs, called_fusions)


#******************************************************************  
# MAP CALLED FUSION BREAKPOINTS TO KNOWN GENES OF .bed FILE
#****************************************************************** 
#known_genes_file_name ="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID_called_genes.bed"
#known_genes_file ="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes.bed"
#known_genes_file_name ="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID.bed"
known_genes_file_name=args.annotation_file
gene_pairs, called_fusions = map_called_fusions_to_known_genes(args.cluster_file, known_genes_file_name) 


new_called_fusions_file = open(args.output_file, 'w+')
## Generate called fusions file contating possible gene pairs
#cluster_type    gene1   gene2  max_split_cnt max_span_cnt sample_type disease tools inferred_fusion_type gene1_on_bnd gene1_close_to_bnd  gene2_on_bnd gene2_close_to_bnd dna_supp  samples chr1 breakpoint_1 chr2 breakpoint_2
new_called_fusions_file.write("#cluster_type	gene1	gene2	max_split_cnt	max_span_cnt	sample_type	disease	tools	inferred_fusion_type	gene1_on_bnd	gene1_close_to_bnd	gene2_on_bnd	gene2_close_to_bnd	dna_supp	samples	chr1	breakpoint_1	chr2	breakpoint_2	gene1_candidates	gene2_candidates	gene1_strands	gene2_strands\n")
for gene_pair,called_fusion in zip(gene_pairs, called_fusions):
    left_mapped_genes = gene_pair[0]
    right_mapped_genes = gene_pair[1]
    called_fusion.gene1_candidates = "|".join([gene[0] for gene in left_mapped_genes]) if len(left_mapped_genes) else "NA"
    called_fusion.gene2_candidates = "|".join([gene[0] for gene in right_mapped_genes]) if len(right_mapped_genes) else "NA"
    called_fusion.gene1_strands += "|".join([gene[2] for gene in left_mapped_genes])
    called_fusion.gene2_strands += "|".join([gene[2] for gene in right_mapped_genes])
    print("|".join([gene[2] for gene in left_mapped_genes]), "|".join([gene[2] for gene in right_mapped_genes]))
    print(called_fusion.gene1_strands, called_fusion.gene2_strands)
    val = called_fusion
    #modified "tools" to include comma-separated list of tools. Assumes fusions have already been merged prior to this script being run
    attribute_list = [called_fusion.cluster_type, called_fusion.gene1, called_fusion.gene2, called_fusion.max_split_cnt, called_fusion.max_span_cnt, called_fusion.sample_type, called_fusion.disease[0], ",".join([tool for tool in called_fusion.tools]), called_fusion.inferred_fusion_type, called_fusion.gene1_on_bnd, called_fusion.gene1_close_to_bnd, val.gene2_on_bnd, val.gene2_close_to_bnd, val.dna_supp, val.samples[0], val.chr1, val.breakpoint_1[0], val.chr2, val.breakpoint_2[0], called_fusion.gene1_candidates, called_fusion.gene2_candidates, called_fusion.gene1_strands, called_fusion.gene2_strands]

    attribute_list = [str(attribute) for attribute in attribute_list]
    new_called_fusions_file.write("\t".join(attribute_list ) + "\n")
new_called_fusions_file.close()
###gene1  gene2   NA      NA      chr1  chr2  start1  end1    start2  end2 


#all_genes = map_called_fusions_to_known_genes(args.called_fusions_file, known_genes_file_name) 
#print( "|".join([gene for gene in all_genes]) )

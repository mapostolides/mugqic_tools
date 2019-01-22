#! /usr/bin/env python

"""
Add the following to validated fusions file: gene names, category of gene, and region of gene within which breakpoints are located
"""

import sys
import argparse
from validation_classes import *
from pygeneann import *
parser = argparse.ArgumentParser()

parser.add_argument('validated_fusions_file', action='store', help='A file containing a set of validated fusions')
parser.add_argument('sample_name', action='store', help='Name of the sample being validated')
 

args = parser.parse_args()


def map_validated_fusions_to_known_genes(validated_fusions_file_name, gene_ann_bed):
    """
    Maps start position of validated fusions to a known gene in the known genes .bed file
    """
    genebed_objects = []
    validated_fusions = []
    for line in open(gene_ann_bed, 'r+'):
        genebed_object = GeneBed(line)
        genebed_objects.append(genebed_object)
    for line in open(validated_fusions_file_name):
        if line.startswith('#'): continue
        fusion = ValidatedFusion(line)
        validated_fusions.append(fusion)

    print("MAPPED GENES FOR VALIDATED FUSIONS")
    gene_pairs = []
    for fusion in validated_fusions:
        # check whole .bed file for candidade gene mappings
        for genebed in genebed_objects:
            if ( ( genebed.start <= fusion.start1 <= genebed.end) or (genebed.start <= fusion.end1 <= genebed.end) ) and (str(fusion.chr1)) == genebed.chr:
                fusion_left_names = [tup[0] for tup in fusion.left]
                if genebed.gene_name not in fusion_left_names:
                    fusion.left.append((genebed.gene_name, genebed.type, genebed.strand, [genebed.start, genebed.end]))
            if ( ( genebed.start <= fusion.start2 <= genebed.end) or (genebed.start <= fusion.end2 <= genebed.end) ) and (str(fusion.chr2)) == genebed.chr:
                fusion_right_names = [tup[0] for tup in fusion.right]
                if genebed.gene_name not in fusion_right_names:
                    fusion.right.append((genebed.gene_name, genebed.type, genebed.strand, [genebed.start, genebed.end]))
        gene_pairs.append((fusion.left, fusion.right))
        print("Fusion gene pairs:", fusion.left, fusion.right)
    return (gene_pairs, validated_fusions)


#******************************************************************  
# MAP VALIDATED FUSION BREAKPOINTS TO KNOWN GENES OF .bed FILE
#****************************************************************** 
#known_genes_file_name ="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID_validated_genes.bed"
#known_genes_file ="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes.bed"
known_genes_file_name ="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes_ensID.bed"
gene_pairs, validated_fusions = map_validated_fusions_to_known_genes(args.validated_fusions_file, known_genes_file_name) 


new_validated_fusions_file = open('/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/stjude_validation/process_validated_fusions_file/' + args.sample_name + 'fusions_valfile-Valid-bed_intersected.tsv', 'w+')
## Generate validated fusions file contating possible gene pairs
for gene_pair,validated_fusion in zip(gene_pairs, validated_fusions):
    left_mapped_genes = gene_pair[0]
    right_mapped_genes = gene_pair[1]
    validated_fusion.gene1_candidates = "|".join([gene[0] for gene in left_mapped_genes])
    validated_fusion.gene2_candidates = "|".join([gene[0] for gene in right_mapped_genes])
    validated_fusion.gene1_strands += "|".join([gene[2] for gene in left_mapped_genes])
    validated_fusion.gene2_strands += "|".join([gene[2] for gene in right_mapped_genes])
    print("|".join([gene[2] for gene in left_mapped_genes]), "|".join([gene[2] for gene in right_mapped_genes]))
    print(validated_fusion.gene1_strands, validated_fusion.gene2_strands)
    val = validated_fusion
    new_validated_fusions_file.write("\t".join([val.gene1, val.gene2, validated_fusion.gene1_strands, validated_fusion.gene2_strands, str(val.chr1), str(val.chr2), str(val.start1), str(val.end1), str(val.start2), str(val.end2), val.sample_name, val.gene1_candidates, val.gene2_candidates ]) + "\n")
new_validated_fusions_file.close()
###gene1  gene2   NA      NA      chr1  chr2  start1  end1    start2  end2 


#all_genes = map_validated_fusions_to_known_genes(args.validated_fusions_file, known_genes_file_name) 
#print( "|".join([gene for gene in all_genes]) )

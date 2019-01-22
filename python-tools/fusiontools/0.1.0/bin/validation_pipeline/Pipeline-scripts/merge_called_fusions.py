import sys
import argparse
from pygeneann import *
import os
import sqlite3
from validation_classes import *


# instantiate parser object
parser = argparse.ArgumentParser()

parser.add_argument('fusion_cluster_file', action='store', help='Good luck! Hope that was helpful')                              
parser.add_argument('sample_name', action='store', help='Name of the sample being validated')
parser.add_argument('output_file', action='store', help='absolute path of the output file')

args = parser.parse_args()

def generate_merged_category_list_by_breakpoints(category_files):
    """
    Generates a category list by merging category fusions based on breakpoints 

    :param category_files: a list of .cluster file paths 
    :param output_file: file being written to, containing merged fusions
    :return: list of merged CategoryFusions objects
    """
    # generate list of merged CategoryFusions objects
    categoryfusion_objects = []
    for category_file in category_files:
        category_file = open(category_file, 'r')
        categoryfusion_objects += [CategoryFusions(line) for line in category_file if not line.startswith('#')]
        print("lenght of categoryfusion_objects",len(categoryfusion_objects))
        category_file.close()
    # create dict of merged fusions with gene pairs as keys
    merged_dict = {}
    #collect duplicate calls(optional)
    #duplicate_dict = {}
    for category_fusion in categoryfusion_objects:
        if (category_fusion.gene1_candidates, category_fusion.gene2_candidates) in merged_dict.keys():
            # compare breakpoints
             if (category_fusion.breakpoint_1[0] == merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)].breakpoint_1[0]  ) and (category_fusion.breakpoint_2[0] == merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)].breakpoint_2[0] ):
                 merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)] = merge_fusions(merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)], category_fusion)
                 # to see deFuse dubplicate calls
                 #duplicate_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)].append(category_fusion)                                                                                    
        else:
            merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)] = category_fusion
            # to see deFuse dubplicate calls
            #duplicate_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)] = [category_fusion]
    return [merged_dict[candidate_gene_pair] for candidate_gene_pair in merged_dict.keys()]

def generate_merged_category_list_by_candidate_genes(category_files):
    """
    Generates a category list by merging category fusions based on breakpoints 

    :param category_files: a list of .cluster file paths 
    :param output_file: file being written to, containing merged fusions
    :return: list of merged CategoryFusions objects
    """
    # generate list of merged CategoryFusions objects
    categoryfusion_objects = []
    for category_file in category_files:
        category_file = open(category_file, 'r')
        categoryfusion_objects += [CategoryFusions(line) for line in category_file if not line.startswith('#')]
        print("lenght of categoryfusion_objects",len(categoryfusion_objects)) 
        category_file.close()
    # create dict of merged fusions with bar-delimited candidate genes as keys
    merged_dict = {}
    num = 0
    for category_fusion in categoryfusion_objects:
        num += 1
        if num%5000==0:
            print("{} fusions processed".format(num))
        #merge fusions based on bar-delimited candidate genes
        if (category_fusion.gene2_candidates, category_fusion.gene1_candidates) in merged_dict:
            #print("inverse match")
            fusion_present_in_dictionary = merged_dict[(category_fusion.gene2_candidates, category_fusion.gene1_candidates)] 
            flipped_fusion = switch_positions(category_fusion)
            #fusion is flipped, so need to reverse order of gene1_candidates and gene2_candidates
            merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)] = merge_fusions(fusion_present_in_dictionary, flipped_fusion) 

        elif (category_fusion.gene1_candidates, category_fusion.gene2_candidates) in merged_dict: 
            #print("forward match")
            fusion_present_in_dictionary = merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)] 
            merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)] = merge_fusions(fusion_present_in_dictionary, category_fusion) 

        else:
            #print("new fusion")
            merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)] = category_fusion 
        #print("fusion added to merged_dict", (category_fusion.gene1_candidates, category_fusion.gene2_candidates))
        #print("merged dict after merge", merged_dict.keys())
    return [merged_dict[candidate_gene_pair] for candidate_gene_pair in merged_dict.keys()]

 #       try:
 #           fusion_present_in_dictionary = merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)]  
 #           merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)] = merge_fusions(fusion_present_in_dictionary, category_fusion)
 #       except:    
 #           fusion_present_in_dictionary = merged_dict[(category_fusion.gene2_candidates, category_fusion.gene1_candidates)] 
 #           
 #       except:     
 #           merged_dict[(category_fusion.gene1_candidates, category_fusion.gene2_candidates)] = category_fusion 
 #


def switch_positions(fusion):
    #switch pos order of: chr, breakpoint, gene_strands, gene, gene_candidates
    fusion.chr1, fusion.chr2 = fusion.chr2, fusion.chr1
    fusion.breakpoint_1, fusion.breakpoint_2 = fusion.breakpoint_2, fusion.breakpoint_1
    fusion.gene1_strands, fusion.gene2_strands = fusion.gene2_strands, fusion.gene1_strands
    fusion.gene1, fusion.gene2 = fusion.gene2, fusion.gene1
    fusion.gene1_candidates, fusion.gene2_candidates = fusion.gene2_candidates, fusion.gene1_candidates
    return fusion

def merge_fusions(fusion1, fusion2):
    """
    Merges the contents of 2 CategoryFusions objects
    """
    fusion1.tools += fusion2.tools
    fusion1.max_split_cnt = max(fusion1.max_split_cnt, fusion2.max_split_cnt)
    fusion1.max_span_cnt = max(fusion1.max_span_cnt, fusion2.max_span_cnt)
    fusion1.breakpoint_1 += fusion2.breakpoint_1
    fusion1.breakpoint_2 += fusion2.breakpoint_2
    return fusion1
    # may need to merge more fields later if more fields of CategoryFusion objects become used.

#**********************************************************************************************************************************
    # write out deFuse call file with duplicates clustered
#    defuse_duplicate_file = "cluster_stats_files/" + args.sample_name + "-defuse_calls_ordered_by_duplicates.cluster"
#    category_list = []
#    for key in duplicate_dict.keys():
#        for fusion in duplicate_dict[key]:
#            category_list.append(fusion)
#    category_file = open(defuse_duplicate_file, "w+")
#    category_file.write("#cluster_type    gene1   gene2  max_split_cnt max_span_cnt sample_type disease tools inferred_fusion_type gene1_on_bnd gene1_close_to_bnd  gene2_on_bnd gene2_close_to_bnd dna_supp  samples chr1 breakpoint_1 chr2 breakpoint_2"  + "\n")
#        
#    for category_fusion in category_list:
#        category_file.write(category_fusion.line + "\n")
#    category_file.close() 
#**********************************************************************************************************************************
 




def generate_merged_total_cluster_file(mode, all_cluster_files):
    """ 
    :param mode: either "gene_names" or "breakpoints"
    """
    # Generate total .cluster file merged based on gene names
    if mode == "candidate_genes":
        merged_total_list = generate_merged_category_list_by_candidate_genes(all_cluster_files)
    elif mode == "breakpoints":
        merged_total_list = generate_merged_category_list_by_breakpoints(all_cluster_files)

    merged_total_cluster_file = open(merged_total_cluster_file_name, 'w+')
    merged_total_cluster_file.write("#cluster_type    gene1   gene2  max_split_cnt max_span_cnt sample_type disease tools inferred_fusion_type gene1_on_bnd gene1_close_to_bnd  gene2_on_bnd gene2_close_to_bnd dna_supp  samples chr1 breakpoint_1 chr2 breakpoint_2"  + "\n")
    for fusion in merged_total_list:
        merged_total_cluster_file.write("\t".join(map(str, [fusion.cluster_type, fusion.gene1, fusion.gene2, fusion.max_split_cnt, fusion.max_span_cnt, fusion.sample_type, ",".join(list(set(fusion.disease))), ",".join(list(set(sorted(fusion.tools)))), fusion.inferred_fusion_type, fusion.gene1_on_bnd, fusion.gene1_close_to_bnd, fusion.gene2_on_bnd, fusion.gene2_close_to_bnd, fusion.dna_supp, ",".join(list(set(fusion.samples))), fusion.chr1, ",".join(list(set([str(bp) for bp in fusion.breakpoint_1]))), fusion.chr2, ",".join(list(set([str(bp) for bp in fusion.breakpoint_2]))), fusion.gene1_candidates, fusion.gene2_candidates, fusion.gene1_strands, fusion.gene2_strands ])) + "\n")
    merged_total_cluster_file.close()


# MAIN

all_cluster_files = [args.fusion_cluster_file]

# Generate total .cluster file merged based on breakpoints 
merged_total_cluster_file_name = args.output_file
#generate_merged_total_cluster_file("breakpoints", all_cluster_files)
generate_merged_total_cluster_file("candidate_genes", all_cluster_files)


# Generate total .cluster file merged based on gene names
#merged_total_cluster_file_name = "cluster_stats_files/merged_gene_names.total.cluster" 
#generate_merged_total_cluster_file("gene_names", all_cluster_files)

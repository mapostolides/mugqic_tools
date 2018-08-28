#! /usr/bin/env python

import sys
import argparse
from pygeneann import *
import os

#import numpy as np
#import pandas as pd

# instantiate parser object
parser = argparse.ArgumentParser()

parser.add_argument('validated_fusions_file', action='store', help='A file containing a set of validated fusions')
parser.add_argument('fusion_cluster_file', action='store', help='Fusion reannation file clustered by head/tail genes. (merged.cff.reann.dnasupp.bwafilter.<seq_len>.cluster file)')

args = parser.parse_args()

     
class ValidatedFusionStats(CategoryFusionStats):
    """
    Takes as input a fusion_cluster_file/category_file, and validated_fusions_file.
    Generates statistics about alse positives, false negatives, and true positives with respect to
    the validated fusions present in the validated_fusions_file.
    """

    def __init__(self, validated_fusions_file, fusion_cluster_file, name):
        CategoryFusionStats.__init__(self, fusion_cluster_file)
        # instantiate instance variables
        self.name = name
        self.cluster_files_dir = "cluster_stats_files"
        # list of true positive CategoryFusion objects
        self.true_positives = [] 
        # lists of gene pair tuples       
        self.validated_fusions = []
        self.false_positives = []
        self.false_negatives = []
        # counts
        self.num_true_positives = 0
        self.num_false_positives = 0
        self.num_false_negatives = 0

        if not os.path.exists(self.cluster_files_dir):
            os.mkdir(self.cluster_files_dir)
        self.__load_validated_fusions_file(validated_fusions_file)

    def __load_validated_fusions_file(self, validated_fusions_file):
        validated_fusions_file = open(validated_fusions_file, "r")
        for line_val in validated_fusions_file:
            gene1_val, gene2_val = line_val.split()[0], line_val.split()[1]
            self.validated_fusions.append( (gene1_val, gene2_val) )
        validated_fusions_file.close()
        print "number of validated fusions: ", len(self.validated_fusions) 

    def find_true_positives(self):
        """
        Intersection of validated_fusions and output_fusions gives true positives
        These are all the validated fusions detected by the pipeline
        """
        # A list of CategoryFusion objects that are true positives detected by pipeline
        true_positives_file_name = os.path.join(self.cluster_files_dir, self.name + ".true_positives.cluster" )
        true_positives_file = open(true_positives_file_name, "w+")
        # create file 
        for validated_gene_pair in self.validated_fusions:
            for category_fusion in self.category_list:
                fusion_gene_pair = (category_fusion.gene1, category_fusion.gene2)
                if (validated_gene_pair == fusion_gene_pair):
                    self.true_positives.append(category_fusion)
                    self.num_true_positives += 1
                    true_positives_file.write(category_fusion.line + "\n")
        true_positives_file.close()
    
    def find_false_positives(self):
        # A list of CategoryFusion objects that are true positives detected by pipeline
        false_positives_file_name = os.path.join( self.cluster_files_dir, self.name + ".false_positives.cluster" )
        false_positives_file = open(false_positives_file_name, "w+")
        # create file
        for category_fusion in self.category_list: 
            fusion_gene_pair = (category_fusion.gene1, category_fusion.gene2)
            if fusion_gene_pair not in self.validated_fusions:
                false_positives_file.write(category_fusion.line + "\n")
                self.num_false_positives += 1
                self.false_positives.append(fusion_gene_pair)
        false_positives_file.close()

    def find_false_negatives(self):
        """
        Set difference of validated_fusions - true_positives gives false negatives 
        Finds gene fusions which are in the provided validated_fusions file, but are not
        detected by the pipeline.
        """
        true_positive_gene_pairs = [(category_fusion.gene1, category_fusion.gene2) for category_fusion in self.true_positives]
        for validated_gene_pair in self.validated_fusions:
            if validated_gene_pair not in true_positive_gene_pairs: 
                self.false_negatives.append(validated_gene_pair) 
                self.num_false_negatives += 1
        # create file of tuples
        false_negatives_file_name = os.path.join( self.cluster_files_dir, self.name + ".false_negatives.txt" ) 
        false_negatives_file = open(false_negatives_file_name, "w+") 
        for gene_pair in self.false_negatives:
            false_negatives_file.write(gene_pair[0] + "\t" + gene_pair[1] + "\n")
            
    def write_tuples_to_file(self, tuple_list):
        """
        Outputs each category fusion object in a human-readable format
        """
        for tup in tuple_list:
            sys.stdout.write(tup)

    def generate_category_file(self, output_file):
        """
        Generates a category file from self.category_list
        """
        category_file = open(output_file, "w+")
        for category_fusion in self.category_list:
            category_file.write(category_fusion.line + "\n")
        
         
    def compare_validated_and_output_fusions(self):
        """
        Compares the fusion pairs in the validated fusion file with the detected fusions in the outputted fusion cluster file
        Outputs results/statistics to an output file
        :param self.validated_fusions: a list of tuples containing fusion gene pairs we would expect to be detected 
                                  by the pipeline with the given input sequence data
        :param self.category_list: file merged.cff.reann.dnasupp.bwafilter.30.cluster, outputted by the pipeline
        """
        self.find_true_positives()
        
        self.find_false_positives()
    
        false_negatives = self.find_false_negatives()

    #test: DIDN"T WORK
    def check_false_positives(self):
        """
        Checks to see if either first or second gene of false positive set matches 
        false negative set
        """
        left_validated_fusion_genes = [gene_pair[0] for gene_pair in self.false_negatives]
        left_false_positive_genes = [gene_pair[0] for gene_pair in self.false_positives]
        #left_match = [gene for gene in left_false_positive_genes if gene in left_validated_fusion_genes ]
        left_match = set(left_false_positive_genes).intersection(set(left_validated_fusion_genes))
        
        right_false_positive_genes = [gene_pair[1] for gene_pair in self.false_positives]
        right_validated_fusion_genes = [gene_pair[1] for gene_pair in self.false_negatives]
        #right_match = [gene for gene in right_false_positive_genes if gene in right_validated_fusion_genes ]
        right_match = set(right_false_positive_genes).intersection(set(right_validated_fusion_genes))

        #[x for x in fusion_list if set(sample_list).intersection(set(x.samples))]
        return left_match, right_match
        

def generate_filtered_category_file(fusion_stats, output_file, tool=None, num=None):
    """
    Generates a category file for a specific tool

    :param fusion_stats: a ValidatedFusionStats object
    :param tool: name of tool we are selecting for
    :return: name of category file generated, same as output_file
    """
    # filter fusions
    if tool: 
        filtered_list = fusion_stats.filter_tools_name(fusion_stats.category_list, tool)
    elif num:
        filtered_list = fusion_stats.filter_tools_num(fusion_stats.category_list, 2)
    # generatfalse_positivese
    #category_file = open(os.path.join(fusion_stats.cluster_files_dir, output_file), "w+")
    category_file = open(output_file, "w+")
    for category_fusion in filtered_list:
        category_file.write(category_fusion.line + "\n")
    category_file.close()

    return output_file

# create ValidatedFusionStats object for overall .cluster file 
total_fusion_stats = ValidatedFusionStats(args.validated_fusions_file, args.fusion_cluster_file, "total")

# Divide cluster file into 4 separate files create new ValidatedFusionStats object for each tool
fusion_stats_objects = []

# generate .cluster files for each tool
tools = ["integrate", "fusionmap", "ericscript"]
for tool in tools:
    #generate category file
    category_file = generate_filtered_category_file(total_fusion_stats, "cluster_stats_files/" + tool + ".cluster", tool=tool)
    #create fusion stats object
    fusion_stats_object = ValidatedFusionStats(args.validated_fusions_file, category_file, tool)
    fusion_stats_objects.append(fusion_stats_object)

# generate .cluster files for 2 or more tools
two_or_more_category_file = generate_filtered_category_file(total_fusion_stats, "cluster_stats_files/" + "two_or_more" + ".cluster", num=2)
two_or_more_fusion_stats = ValidatedFusionStats(args.validated_fusions_file, two_or_more_category_file, "two_or_more")
fusion_stats_objects.append(two_or_more_fusion_stats)

# append total last
fusion_stats_objects.append(total_fusion_stats)

# generate stats
print 'Tool	num_fusions	num_true_pos	num_false_pos	num_false_neg'
for fusion_stats in fusion_stats_objects:
    #print "generating cluster files for " + fusion_stats_object.name
    # generate TP, FP and FN files
    fusion_stats.compare_validated_and_output_fusions()
    print '{}	{}	{}	{}	{}	'.format(fusion_stats.name, fusion_stats.num_fusions, fusion_stats.num_true_positives, 
                                                         fusion_stats.num_false_positives, fusion_stats.num_false_negatives)
#DIDN"T WORK
# check for differently named genes: check to see if at least one of the false_neg genes matches false negative set
integrate_left, integrate_right = fusion_stats_objects[0].check_false_positives()

print "left false +ives that match left false -ives: ", integrate_left
print "right false +ives that match right false -ives: ",integrate_right


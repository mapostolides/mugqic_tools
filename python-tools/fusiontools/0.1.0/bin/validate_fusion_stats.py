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

class ValidatedFusion():
    """
    Represents one line of a validated fusion. Takes in a line in the following format:
    gene1   gene2   ensg1   ensg2   chr1    chr2    start1  end1    start2  end2
    """
    def __init__(self, validated_fusion_line):
        self.__load_validated_fusion(validated_fusion_line)
    
    def __load_validated_fusion(self, validated_fusion_line): 
        line = validated_fusion_line.split()
        self.gene1 = line[0]
        self.gene2 =line[1]
        self.ensg1 = line[2]
        self.ensg2 = line[3]
        self.chr1 = line[4]
        self.chr2 = line[5]
        self.start1 = int(line[6])
        self.end1 = int(line[7])
        self.start2 = int(line[8])
        self.end2 = int(line[9])
        self.line = validated_fusion_line
     
class ValidatedFusionStats(CategoryFusionStats):
    """
    Takes as input a fusion_cluster_file/category_file, and validated_fusions_file.
    Generates statistics about alse positives, false negatives, and true positives with respect to
    the validated fusions present in the validated_fusions_file.
    """

    def __init__(self, validated_fusions_file, fusion_cluster_file, name):
        # cluster file is loaded into self.category_list in CategoryFusionStats class definition
        CategoryFusionStats.__init__(self, fusion_cluster_file)
        # instantiate instance variables
        self.name = name
        self.cluster_files_dir = "cluster_stats_files"
        # list of true positive CategoryFusion objects
        self.true_positives = []
        # list of ValidatedFusion objects
        self.validated_fusions = []
        self.true_positive_validated_fusions = []
        self.false_negative_validated_fusions = []
        # lists of gene pair tuples       
        self.false_positives = []
        self.true_positive_gene_pairs = []
        self.validated_gene_pairs = [] 
        # counts
        self.num_true_positives = 0
        self.num_unvalidated_fusions = 0 
        self.num_false_negatives = 0 
        #make cluster_stats_files directory if it doesn't already exist
        if not os.path.exists(self.cluster_files_dir):
            os.mkdir(self.cluster_files_dir)
        self.__load_validated_fusions_file(validated_fusions_file)
        
        # load category file, using altrenate __load

    def __load_validated_fusions_file(self, validated_fusions_file):
        validated_fusions_file = open(validated_fusions_file, "r")
        for line in validated_fusions_file:
            if line.startswith('#'): continue
            fusion = ValidatedFusion(line)
            self.validated_fusions.append(fusion)        
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
        # unvalidated fusions
        self.unvalidated_fusions = self.category_list
        for fusion in self.validated_fusions:
            validated_gene_pair = (fusion.gene1, fusion.gene2)
            for category_fusion in self.category_list:
                # check if fusion breakpoints are within genomic loci
                if ( (fusion.start1 <= min(category_fusion.breakpoint_1) and min(category_fusion.breakpoint_1) <= fusion.end1) and 
                    ( fusion.start2 <= max(category_fusion.breakpoint_2) and max(category_fusion.breakpoint_2) <= fusion.end2)): 

                    self.true_positive_gene_pairs.append(validated_gene_pair)
                    self.true_positives.append(category_fusion)
                    # remove fusions that are true positives to generate unvalidated_fusions list
                    self.unvalidated_fusions.remove(category_fusion)
                    if fusion not in self.true_positive_validated_fusions:
                        self.true_positive_validated_fusions.append(fusion)
                    true_positives_file.write(category_fusion.line + "\n")              
        true_positives_file.close()
        self.num_true_positives = len(self.true_positive_gene_pairs)

    #def compare_gene_names():
    #    compare_gene_names_file_name = os.path.join(self.cluster_files_dir, self.name + ".compare_gene_names_file.txt" )
        #compare_gene_names_file = open(compare_gene_names_file_name, "w+")
        #compare_gene_names_file.write("fusion.gene1    fusion.gene2    category_fusion.gene1   category_fusion.gene2\n")
                    #to_write = "{}	{}	{}	{}\n".format(fusion.gene1, fusion.gene2, category_fusion.gene1, category_fusion.gene2) 
                    #compare_gene_names_file.write(to_write)
        #compare_gene_names_file.close()

    def find_unvalidated_fusions(self):
        # A list of CategoryFusion objects that are unvalidated detected by pipeline
        # False positives and possible true fusions that have not been validated
        unvalidated_fusions_file_name= os.path.join( self.cluster_files_dir, self.name + ".unvalidated_fusions.cluster" )
        unvalidated_fusions_file = open(unvalidated_fusions_file_name, "w+")
        # create file
        #self.unvalidated_fusions = [fusion for fusion in self.category_list if fusion not in [true_pos for true_pos in self.true_positives]]
        for fusion in self.unvalidated_fusions: 
                unvalidated_fusions_file.write(fusion.line + "\n")
        unvalidated_fusions_file.close()
        self.num_unvalidated_fusions = len(self.unvalidated_fusions) 

    def find_false_negatives(self):
        """
        Set difference of validated_fusions - true_positives gives false negatives 
        Finds gene fusions which are in the provided validated_fusions file, but are not
        detected by the pipeline.
        """
        self.false_negative_validated_fusions = [validated_fusion for validated_fusion in self.validated_fusions if validated_fusion not in self.true_positive_validated_fusions]  
        false_negatives_file_name = os.path.join( self.cluster_files_dir, self.name + ".false_negatives.txt" ) 
        false_negatives_file = open(false_negatives_file_name, "w+") 
        for fusion in self.false_negative_validated_fusions:
            false_negatives_file.write(fusion.line)
        false_negatives_file.close()
        self.num_false_negatives = len(self.false_negative_validated_fusions)
 
    def calculate_sensitivity(self):
        self.sensitivity = float(self.num_true_positives)/float(len(self.validated_fusions))

    def calculate_precision(self):
        self.precision = float(self.num_true_positives)/float(self.num_true_positives + self.num_unvalidated_fusions)

    def compare_validated_and_output_fusions(self):
        """
        Compares the fusion pairs in the validated fusion file with the detected fusions in the outputted fusion cluster file
        Outputs results/statistics to an output file
        :param self.validated_fusions: a list of tuples containing fusion gene pairs we would expect to be detected 
                                  by the pipeline with the given input sequence data
        :param self.category_list: file merged.cff.reann.dnasupp.bwafilter.30.cluster, outputted by the pipeline
        """
        self.find_true_positives()
        
        self.find_unvalidated_fusions()
    
        self.find_false_negatives()
        
        self.calculate_sensitivity()

        self.calculate_precision()

    def convert_cff_to_fake_category_file(cff_file):
        """
        Converts cff file into category file format, inserting dummy values where needed so that file can be processed
        by validation tool
        """
        #command =         

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
        filtered_list = fusion_stats.filter_tools_num(fusion_stats.category_list, num)
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
tools = ["defuse", "integrate", "fusionmap", "ericscript"]
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
print 'Tool	num_fusions	num_true_pos	num_unvalidated	num_false_neg	sensitivity	precision'
for fusion_stats in fusion_stats_objects:
    #print "generating cluster files for " + fusion_stats_object.name
    # generate TP, FP and FN files
    fusion_stats.compare_validated_and_output_fusions()
    #print '{}	{}	{}'.format(fusion_stats.name, fusion_stats.num_fusions, fusion_stats.num_true_positives) 
    print '{}	{}	{}	{}	{}	{}	{}	'.format(fusion_stats.name, fusion_stats.num_fusions, fusion_stats.num_true_positives, 
                                                         fusion_stats.num_unvalidated_fusions, fusion_stats.num_false_negatives, fusion_stats.sensitivity, fusion_stats.precision)
 
# looks at true positives for total_fusions. Some are repeated, but duplication is not detected by the pipeline
for category_fusion in total_fusion_stats.true_positives:
    print(category_fusion.tools, category_fusion.chr1, category_fusion.chr2, category_fusion.breakpoint_1, category_fusion.breakpoint_2)

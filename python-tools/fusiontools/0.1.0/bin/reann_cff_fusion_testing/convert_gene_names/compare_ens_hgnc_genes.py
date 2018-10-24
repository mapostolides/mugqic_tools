#!/usr/bin/env python

import subprocess
import pandas as pd
import os
testing = 1

# from: gene_fusion reference .bed file obtained by Yue somehow
ens_gene_file = open("/Users/mapostolides/mugqic_tools-debugging/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/ens_known_genes.bed")

# obtained from [https://www.genenames.org/cgi-bin/download]
#hgnc_gene_file = open("hgnc_complete_set.txt")
hgnc_gene_file = open("/Users/mapostolides/mugqic_tools-debugging/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/convert_gene_names/gene_with_protein_product.txt")

# from: [ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz]
# tab-delimited file with a column contating all aliases for each gene
ncbi_gene_info_file = open("/Users/mapostolides/mugqic_tools-debugging/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/convert_gene_names/Homo_sapiens.gene_info") 

# file containing candidate fusions from all 4 fusion callers 
sim45_fusion_output = open("/Users/mapostolides/mugqic_tools-debugging/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/convert_gene_names/smc_rna_sim45_outputOct5_merged-truncated.cff")

# generate left and right gene lists from original .cff file
table = pd.read_table(sim45_fusion_output, header=None)
left_genes = [item[0] for item in table.iloc[0:, 13:14].values.tolist()]
#left_genes = list(table.iloc[0:, 13:14])#.tolist()
right_genes = [item[0] for item in table.iloc[0:, 15:16].values.tolist()]
#right_genes = table.iloc[0:, 15:16]

# convert left and right genes from original .cff file to "limma" R package gene names
# run R script here
arg = '/Users/mapostolides/mugqic_tools-debugging/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/convert_gene_names/smc_rna_sim45_outputOct5_merged-truncated.cff'

# Run external R script, store output using pipe
print("RUNNING R limma SCRIPT")
p = subprocess.Popen('RScript convert_genes_limma.R ' + arg, stdout=subprocess.PIPE,  shell=True)
(output, err) = p.communicate()
print("R SUBPROCESS COMPLETE")

# format R output
output = str(output).split("\"")
left_genes_renamed, right_genes_renamed = output[1].split(), output[3].split()



# replace NA values 
for i in range(0,len(left_genes_renamed)):
    if left_genes_renamed[i] == 'NA':
        left_genes_renamed[i] = right_genes[i]
    if right_genes_renamed[i] == 'NA':
        right_genes_renamed[i] = right_genes[i]
print("before: ", left_genes)
print("after: ", left_genes_renamed)
print("before: ", right_genes)
print("after: ", right_genes_renamed)

renamed_gene_file = open(arg + ".renamed", 'w+')
#  write new file
i=0
for line in open(arg):
    line=line.split()
    line[13], line[15] = left_genes_renamed[i], right_genes_renamed[i]
    line = "\t".join(str(x) for x in line)
    print(line)
    renamed_gene_file.write(line+"\n")
    i+=1
renamed_gene_file.close()


# SCRAP BELOW

#print("orig_len:", len(left_genes), "renamed_len:", len(left_genes_renamed) )
#print(left_genes_renamed, right_genes_renamed)

#p = subprocess.check_output(['RScript convert_genes_limma.R', arg],  shell=True)




# make dictionary of {consensus_name : [aliases] } using __ file
if testing == 2:
    intersection = ens_gene_set.intersection(hgnc_gene_set)
    #print(ens_gene_set)
    print("ens gene set length:", len(ens_gene_set))
    print("hgnc gene set length:", len(hgnc_gene_set))
    print("intersection length:", len(intersection))
    print("ens_gene_set - hgnc_gene_set:", ens_gene_set.difference(hgnc_gene_set))
    # get set of genes in ens_gene_file
    ens_genes= []
    for line in ens_gene_file:
        line = line.split()
        gene_name = line[7]
        ens_genes.append(gene_name)
    ens_gene_file.close()
    ens_gene_set = set(ens_genes)
    
    # get set of genes in hgnc_gene_file
    hgnc_genes = []
    for line in hgnc_gene_file:
        line = line.split()
        gene_name = line[1]
        hgnc_genes.append(gene_name)
    hgnc_gene_file.close()
    hgnc_gene_set = set(hgnc_genes)
    

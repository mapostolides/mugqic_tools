#!/usr/bin/env python
#/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts/rename_cff_file_genes-GENAP.py
import os
import sys
import  pygeneann_reads_capture_DEV as pygeneann
import pandas as pd
import sequtils
import pysam
import argparse



#parser = argparse.ArgumentParser()

#parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')



#args = parser.parse_args()


# tab-delimited file with a column contating all aliases for each gene
#tax_id	GeneID	Symbol	LocusTag	Synonyms	dbXrefs	chromosome	map_location	description	type_of_gene	Symbol_from_nomenclature_authority	Full_name_from_nomenclature_authority	Nomenclature_status	Other_designations	Modification_date	Feature_type
#9606	1	A1BG	-	A1B|ABG|GAB|HYST2477	MIM:138670|HGNC:HGNC:5|Ensembl:ENSG00000121410	19	19q13.43	alpha-1-B glycoprotein	protein-coding	A1BG	alpha-1-B glycoprotein	O	alpha-1B-glycoprotein|HEL-S-163pA|epididymis secretory sperm binding protein Li 163pA	20200313	-

#Open NCBI file and create df
ncbi_gene_info_file = open("/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/Homo_sapiens.gene_info") 
df = pd.read_csv(ncbi_gene_info_file, sep='\t')
#df = df.head(4)
#df = df.tail(4)

#https://stackoverflow.com/questions/26336251/pandas-rename-single-dataframe-column-without-knowing-column-name
#Create one row for each alias
df = pd.concat([pd.Series(row['Symbol'], row['Synonyms'].split('|'))
	for _, row in df.iterrows()]).reset_index()
df.rename(columns = {list(df)[0]: 'Alias', list(df)[1]: 'HGNC_Symbol'}, inplace = True)

def alias2hgnc(df, query):
    #check to see if query is already HGNC symbol
    hgnc = df.loc[df.HGNC_Symbol == query].HGNC_Symbol.values.tolist()
    is_hgnc = True if len(hgnc) > 0 else False 
    if is_hgnc: 
        return query
    else:
        hgnc_lst = df.loc[df.Alias == query].HGNC_Symbol.values.tolist()
        return "NA" if len(hgnc_lst) == 0 else hgnc_lst[0]
    
#cff_file = args.cff_file
cff_file ="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/reann_cff_fusion_testing/FOSB--AADACL2/FOSB--AADACL2.cff.MOD"
for line in open(cff_file, "r"):
    fusion = pygeneann.CffFusion(line)
    #handle list of genes from caller
    gene1_lst = fusion.t_gene1.split(",")
    gene1_hgnc = []
    for gene in gene1_lst:
        hgnc_gene = alias2hgnc(df, fusion.t_gene1)
        fusion.t_gene1 = hgnc_gene if hgnc_gene != "NA" else fusion.t_gene1
    
    print(fusion.tostring())    
    #gene2_lst = fusion.t_gene2.split(",")
    #hgnc_gene1 = alias2hgnc(df, fusion.t_gene1) 
    #hgnc_gene2 = alias2hgnc(df, fusion.t_gene2) 
    #print( hgnc_gene1,hgnc_gene2)

## replace NA values with original values to prevent lost information 
#for i in range(0,len(left_genes_renamed)):
#    if left_genes_renamed[i] == 'NA':
#        #left_genes_renamed[i] = 'NA' 
#        left_genes_renamed[i] = left_genes[i]
#    if right_genes_renamed[i] == 'NA': 
#        #right_genes_renamed[i] = 'NA' 
#        right_genes_renamed[i] = right_genes[i]
#
##sys.stderr.write("LEFT GENES RENAMED REPLACED NA"+str(left_genes_renamed) + "\n")
#
#renamed_gene_file = open(cff_file+ ".renamed", 'w+')
##  write new file
#i=0
#for line in open(cff_file):
#    line=line.split()
#    line[13], line[15] = left_genes_renamed[i], right_genes_renamed[i]
#    line = "\t".join(str(x) for x in line)
#    # print line for output to modified .cff file
#    print(line)
#    renamed_gene_file.write(line+"\n")
#    i+=1
#renamed_gene_file.close()
#

# SCRAP BELOW

#print("orig_len:", len(left_genes), "renamed_len:", len(left_genes_renamed) )
#print(left_genes_renamed, right_genes_renamed)

#p = subprocess.check_output(['RScript convert_genes_limma.R', arg],  shell=True)



testing=0
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
    

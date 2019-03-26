import sys
import argparse
from pygeneann import *
import os
from validation_classes import *
import pandas as pd

# INSTANTIATE ARGPARSER AND DEFINE ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument('FI_output_file', action='store', help='FusionInspector output file')
parser.add_argument('renamed_cff_file', action='store', help='.cff file with genes renamed using "Limma" ')
parser.add_argument('FI_validated_cff_file', action='store', help='to write, .cff file with FusionInspector-validated cff fusions')
parser.add_argument('FI_unvalidated_cff_file', action='store', help='to write, .cff file with cff fusions NOT detected by FusionInspector')
parser.add_argument('FI_summary_file', action='store', help='Fusioninspector summary file for all fusion callers ')
parser.add_argument('caller_name', action='store', help='name of caller')

#parser.add_argument('output_file', action='store', help='absolute path of the output file')
args = parser.parse_args()

# GENERATE LIST OF FUSIONINSPECTOR-VALIDATED GENE FUSION PAIRS
FI_output_file = open(args.FI_output_file, 'r+')
FI_output_fusions = []
FI_output_genepair_bps = []
for line in FI_output_file:
    if line.split()[0] == "#FusionName":
        continue
    FI_fusion = FusionInspectorFusion(line)
    FI_output_fusions.append(FI_fusion)
    FI_output_genepair_bps.append((FI_fusion.gene_pair, FI_fusion.bp1, FI_fusion.bp2))
FI_output_file.close()

def is_same_fusion(cff_fusion, FI_fusion):
    """
    Checks whether breakpoints AND gene names match for the two fusions
    """
    if (FI_fusion.bp1 == cff_fusion.pos1 and FI_fusion.bp2 == cff_fusion.pos2) and (FI_fusion.gene1 == cff_fusion.t_gene1 and FI_fusion.gene2 == cff_fusion.t_gene2):
    #if (FI_fusion.gene1 == cff_fusion.t_gene1 and FI_fusion.gene2 == cff_fusion.t_gene2):
        return True
    return False


# SORT CFF FUSIONS INTO FI-VALIDATED AND FI-UNVALIDATED 
renamed_cff_file = open(args.renamed_cff_file, 'r+')
FI_validated_cff_fusions = []
all_cff_fusions = []
for line in renamed_cff_file:
    cff_fusion = CffFusion(line)
    #gene_fusion = cff_fusion.t_gene1 + "--" + cff_fusion.t_gene2
    #all_gene_fusions.append(gene_fusion)
    all_cff_fusions.append(cff_fusion)
    #GATHER VALIDATED FI FUSIONS
    for FI_fusion in FI_output_fusions:
        if is_same_fusion(cff_fusion, FI_fusion): 
            FI_validated_cff_fusions.append(cff_fusion) 
    #REMOVE UNVALIDATED FUSIONS FROM all_cff_fusions
    FI_unvalidated_cff_fusions = set(all_cff_fusions).difference(set(FI_validated_cff_fusions))
renamed_cff_file.close()

#WRITE FI-VALIDATED CFF FUSIONS TO FILE
FI_validated_cff_file = open(args.FI_validated_cff_file, 'w+')
#  write new file
for cff_fusion in FI_validated_cff_fusions:
    FI_validated_cff_file.write(cff_fusion.line)
FI_validated_cff_file.close()

#WRITE FI-UNVALIDATED CFF FUSIONS TO FILE
FI_unvalidated_cff_file = open(args.FI_unvalidated_cff_file, 'w+')
#  write new file
for cff_fusion in FI_unvalidated_cff_fusions:
    FI_unvalidated_cff_file.write(cff_fusion.line)
FI_unvalidated_cff_file.close()

#WRITE SUMMARY STATS TO SUMMARY FILE
FI_summary_file = open(args.FI_summary_file, 'a')
FI_summary_file.write(args.caller_name + "\t" + str(len(all_cff_fusions)) + "\t" + str(len(FI_output_fusions)) + "\t" + str(len(FI_validated_cff_fusions)) + "\t" + str(len(FI_unvalidated_cff_fusions)) + "\n")
FI_summary_file.close()

testing = 0
if testing:

    print("FI OUTPUT FUSIONS")
    print([fusion.gene1 + "--" + fusion.gene2 for fusion in FI_output_fusions])
    #print(FI_output_genepair_bps)
    print("NUMBER OF FI OUTPUT FUSIONS", len(FI_output_fusions))
    
    print("FUSIONS VALIDATED BY FI")
    #FI_validated_gene_pairs=set([fusion.t_gene1 + "--" + fusion.t_gene2 for fusion in FI_validated_cff_fusions])
    FI_validated_gene_pairs=[fusion.t_gene1 + "--" + fusion.t_gene2 for fusion in FI_validated_cff_fusions]
    print(FI_validated_gene_pairs)
    print("NUMBER OF VALIDATED FI FUSIONS", len(FI_validated_gene_pairs))
    
    print("UNVALIDATED CALLER CALLS:", [fusion.t_gene1 + "--" + fusion.t_gene2 for fusion in FI_unvalidated_cff_fusions])
    print("NUMBER OF UNVALIDATED FUSIONS:", len(FI_unvalidated_cff_fusions))

# SCRAP

# generate left and right gene lists from original .cff file
#cff_table = pd.read_table(renamed_cff_file, header=None)
#cff_left_genes = [item[0] for item in cff_table.iloc[0:, 13:14].values.tolist()]
#left_genes = list(table.iloc[0:, 13:14])#.tolist()
#cff_right_genes = [item[0] for item in cff_table.iloc[0:, 15:16].values.tolist()]




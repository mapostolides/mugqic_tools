import os

# INSTANTIATE ARGPARSER AND DEFINE ARGUMENTS
#parser = argparse.ArgumentParser()
#parser.add_argument('cluster_files_dir', action='store', help='Directory containing cluster files, and false negative files')

#parser.add_argument('renamed_cff_file', action='store', help='.cff file with genes renamed using "Limma" ')
#parser.add_argument('FI_validated_cff_file', action='store', help='to write, .cff file with FusionInspector-validated cff fusions')

#args = parser.parse_args()

#cluster_files_dir = args.cluster_files_dir
#cluster_files_dir = "/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/sim45_validation/output-2019-02-04/cluster_stats_files"
cluster_files_dir = "/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/sim52_validation/output-2019-01-25/cluster_stats_files"
#cluster_files_dir = "/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/stjude_validation/output-2019-01-24/cluster_stats_files"

tool_list = ["star_fusion", "defuse", "fusionmap", "ericscript", "integrate"]

files = os.listdir(cluster_files_dir)
single_caller_false_neg_files = [ fil for fil in files if "false_negatives" in fil and "two_or_more" not in fil and "swp" not in fil and any(tool in fil for tool in tool_list) ]
print(single_caller_false_neg_files) 

#CREATE DICTIONARY OF VALIDATED FUSIONS: { genebed_intersection_list : num_callers_missed }
false_neg_dict = {}
for fil in single_caller_false_neg_files:
    file_path = fil
    fil = os.path.join(cluster_files_dir, fil)
    print(fil)
    fil = open(fil, 'r+')
    for line in fil:
        line = line.split()
        if line[0] == "gene1": continue
        if (line[11], line[12]) not in false_neg_dict:
            false_neg_dict[(line[11], line[12])] = [1, [file_path.split(".")[0]]]
        else:
            false_neg_dict[(line[11], line[12])][0] += 1
            false_neg_dict[(line[11], line[12])][1].append(file_path.split(".")[0])
    fil.close()
num = 2 
for tup in [(gene_pair, false_neg_dict[gene_pair]) for gene_pair in false_neg_dict.keys() if false_neg_dict[gene_pair] >= num]:
#for tup in [(gene_pair, false_neg_dict[gene_pair]) for gene_pair in false_neg_dict.keys() if false_neg_dict[gene_pair] == num]:
    print(str(tup[0][0]) + "--" + str(tup[0][1]) +  ":" + str(tup[1]))
        

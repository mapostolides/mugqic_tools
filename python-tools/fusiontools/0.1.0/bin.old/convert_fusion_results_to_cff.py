#!/usr/bin/env python

import sys
import os
import argparse
class FusionResult():
	def __init__(self, tool, line, idxes):
		tmp = line.split()
		idx_chr1, idx_pos1, idx_strand1, idx_chr2, idx_pos2, idx_strand2, idx_split_cnt, idx_pair_cnt, idx_gene1, idx_gene2, idx_gene_location1, idx_gene_location2 = idxes

		self.tool = tool
		self.chr1 = tmp[idx_chr1]
		self.pos1 = tmp[idx_pos1]
		self.strand1 = tmp[idx_strand1] if idx_strand1 != "NA" else "NA" 
		self.chr2 = tmp[idx_chr2]
		self.pos2 = tmp[idx_pos2]
		self.strand2 = tmp[idx_strand2] if idx_strand2 != "NA" else "NA"
		self.split_cnt = tmp[idx_split_cnt] if idx_split_cnt != "NA" else -1
		self.pair_cnt = tmp[idx_pair_cnt] if idx_pair_cnt != "NA" else -1
		self.gene1 = tmp[idx_gene1]
		self.gene2 = tmp[idx_gene2]
		self.gene_location1 = tmp[idx_gene_location1] if idx_gene_location1 != "NA" else "NA"
		self.gene_location2 = tmp[idx_gene_location2] if idx_gene_location2 != "NA" else "NA"
class FusionResultFile():
	def __init__(self, result_file):
		self.fusion_results = []
		for line in open(result_file, "r"):
			tmp = line.split()
			if tmp[0] == "cluster_id": ## defuse header
				self.tool = "Defuse"
			
				self._idx_chr1 = tmp.index("gene_chromosome1")
				self._idx_chr2 = tmp.index("gene_chromosome2")
				self._idx_pos1 = tmp.index("genomic_break_pos1")
				self._idx_pos2 = tmp.index("genomic_break_pos2")
				self._idx_strand1 = tmp.index("genomic_strand1")
				self._idx_strand2 = tmp.index("genomic_strand2")
				self._idx_split_cnt = tmp.index("splitr_count")
				self._idx_pair_cnt = tmp.index("span_count")
				self._idx_gene1 = tmp.index("gene_name1")
				self._idx_gene2 = tmp.index("gene_name2")
				self._idx_gene_location1 = tmp.index("gene_location1")
				self._idx_gene_location2 = tmp.index("gene_location2")
			elif tmp[0] == "FusionID": ## fusionmap header
				self.tool = "FusionMap"

				self._idx_chr1 = tmp.index("Chromosome1")
				self._idx_chr2 = tmp.index("Chromosome2")
				self._idx_pos1 = tmp.index("Position1")
				self._idx_pos2 = tmp.index("Position2")
				self._idx_strand1 = "NA"
				self._idx_strand2 = "NA"
				for i in range(len(tmp)):
					name = tmp[i]
					if "UniqueCuttingPositionCount" in name:
						self._idx_split_cnt = i
						break
				self._idx_pair_cnt = "NA"
				self._idx_gene1 = tmp.index("KnownGene1")
				self._idx_gene2 = tmp.index("KnownGene2")
				self._idx_gene_location1 = "NA"
				self._idx_gene_location2 = "NA"
			elif tmp[0] == "GeneName1": ## ericscript header
				self.tool = "EricScript"

				self._idx_chr1 = tmp.index("chr1")
				self._idx_chr2 = tmp.index("chr2")
				self._idx_pos1 = tmp.index("Breakpoint1")
				self._idx_pos2 = tmp.index("Breakpoint2")
				self._idx_strand1 = "NA"
				self._idx_strand2 = "NA"
				self._idx_split_cnt = tmp.index("spanningreads")
				self._idx_pair_cnt = tmp.index("crossingreads")
				self._idx_gene1 = tmp.index("GeneName1")
				self._idx_gene2 = tmp.index("GeneName2")
				self._idx_gene_location1 = "NA"
				self._idx_gene_location2 = "NA"

			elif tmp[0] == "#5P": ## integrate header
				self.tool = "Integrate"

				self._idx_chr1 = tmp.index("Chr1")
				self._idx_chr2 = tmp.index("Chr2")
				self._idx_pos1 = tmp.index("RNA_BK1")
				self._idx_pos2 = tmp.index("RNA_BK2")
				self._idx_strand1 = "NA"
				self._idx_strand2 = "NA"
				self._idx_split_cnt = "NA"
				self._idx_pair_cnt = "NA"
				self._idx_gene1 = tmp.index("#5P")
				self._idx_gene2 = tmp.index("3P")
				self._idx_gene_location1 = "NA"
				self._idx_gene_location2 = "NA"
					
			else:
				fusion = FusionResult(self.tool, line, [self._idx_chr1, self._idx_pos1, self._idx_strand1, self._idx_chr2, self._idx_pos2, self._idx_strand2, self._idx_split_cnt, self._idx_pair_cnt, self._idx_gene1, self._idx_gene2, self._idx_gene_location1, self._idx_gene_location2])
				if fusion.pos1.isdigit() and fusion.pos2.isdigit(): 
					self.fusion_results.append(fusion)

class SampleInfo():
	def __init__(self, line, idxes):
		tmp = line.split()
		idx_sample, idx_disease, idx_lib, idx_sample_type = idxes
		self.args.sample = tmp[idx_sample]
		self.disease = tmp[idx_disease]
		self.lib = tmp[idx_lib]
		self.sample_type = tmp[idx_sample_type]
class SampleInfoFile():
	def __init__(self, sampleinfo_file):
		self.sampleinfo_dict ={}
		for line in open(sampleinfo_file, "r"):
			tmp = line.split()
			if tmp[0] == "Sample": ## header
				self._idx_sample = tmp.index("Sample")
				self._idx_disease = tmp.index("Disease")
				self._idx_lib = tmp.index("Library")
				self._idx_sample_type = tmp.index("Sample_type")
			else:
				idxes = [self._idx_sample, self._idx_disease, self._idx_lib, self._idx_sample_type]
				self.sampleinfo_dict.setdefault(args.sample, SampleInfo(line, idxes))



parser = argparse.ArgumentParser()
parser.add_argument("sample", help="Sample name")
parser.add_argument("sample_type", help="Tumor or Noraml ", choices=["Tumor", "Normal"])
parser.add_argument("disease_name", help="Disease name")
parser.add_argument("tool", help="Tool name")
parser.add_argument("fusion_result_file", help="Original fusion result file generated by tools")
parser.add_argument("out_dir", help="Output folder")
args = parser.parse_args()

'''				
if len(sys.argv) != 7:
	print "convert_fusion_results_to_cff.py sample sample_type disease_name tool fusion_results outdir"
	sys.exit(1)
sample = sys.argv[1]
sample_type = sys.argv[2]
disease_name = sys.argv[3]
tool = sys.argv[4]
fusion_result_file = sys.argv[5]
out_dir = sys.argv[6]
'''
fusion_results = FusionResultFile(args.fusion_result_file).fusion_results
'''
## CFF format
#Breakpoints Zone
# if strand info missing (NA), it is assumed that gene1 is 5 prime gene and gene2 is 3 prime gene, strand will be infered based on this assumption.

1 chromosome1
2 break_pos1
3 strand1
4 chromosome2
5 break_pos2
6 strand2

## Sample Info Zone

7 data_type  (DNA/RNA)
8 sample_id
9 sample_type (Tumor/Normal)
10 disease

## Software Zone
#
11 tool
12 split_count
13 spanning_count

## Original Annotation Zone

14 gene1
15 area_type1 (cds/intron/utr)
16 gene2
17 area_type2 (cds/intron/utr)
'''
out_file = open(os.path.join(args.out_dir, args.sample + "." + args.tool + ".cff"), "w")
for fusion in fusion_results:
	
	print >> out_file, "\t".join(map(str, [fusion.chr1, fusion.pos1, fusion.strand1, fusion.chr2, fusion.pos2, fusion.strand2, "RNA", args.sample, args.sample_type, args.disease_name, args.tool, fusion.split_cnt, fusion.pair_cnt, \
			fusion.gene1, fusion.gene_location1, fusion.gene2, fusion.gene_location2]))

out_file.close()
		






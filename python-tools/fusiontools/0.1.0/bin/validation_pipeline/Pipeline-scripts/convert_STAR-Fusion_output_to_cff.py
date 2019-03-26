import argparse

#FusionName     JunctionReadCount       SpanningFragCount       SpliceType      LeftGene        LeftBreakpoint  RightGene       RightBreakpoint LargeAnchorSupport      FFPM    LeftBreakDinuc  LeftBreakEntropy        RightBreakDinuc RightBreakEntropy       annots
#PSMC5--ENY2     1080    960     ONLY_REF_SPLICE PSMC5^ENSG00000087191.8 chr17:61906923:+        ENY2^ENSG00000120533.8  chr8:110351549:+        YES_LDAS        34      GT      1.7232  AG      1.9899  ["INTERCHROMOSOMAL[chr17--chr8]"]

# Create STAR_Fusion dict

# instantiate parser object
parser = argparse.ArgumentParser()

parser.add_argument('star_fusion_outfile', action='store', help='A file containing candidate fusions called by STAR-Fusion')
parser.add_argument('star_fusion_cff', action='store', help='A file being written to, in .cff format, containing candidate fusions called by STAR-Fusion')

args = parser.parse_args()

#STAR_Fusion_outfile_path="/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/STAR-Fusion/STAR-Fusion_outputDec6-v1.5.0/star-fusion.fusion_predictions.abridged.tsv"
#STAR_Fusion_outfile_path="/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/STAR-Fusion/STAR-Fusion_outputDec13-v1.5.0-skewer_trimmed_reads/star-fusion.fusion_predictions.abridged.tsv"
STAR_Fusion_outfile_path= args.star_fusion_outfile
STAR_Fusion_outfile = open(STAR_Fusion_outfile_path, 'r+')

lines=[]
for line in STAR_Fusion_outfile:
   lines.append(line)
STAR_Fusion_outfile.close()

cff_file_path = args.star_fusion_cff 
cff_file = open(cff_file_path, 'w+')
for line in lines:
    tmp = line.split()
    if "#FusionName" in line: continue

    STAR_Fusion_outfile.close()# Create STAR_Fusion dict
    #19      14779663        +       19      14800720        -       RNA     smc_rna_sim45   Tumor   VALIDATION      defuse  2       4       EMR3    coding  ZNF333  utr5p
    tool = "star_fusion"
    left_chr_bp_strand = tmp[5].split(":") #chr17:61906923:+
    right_chr_bp_strand = tmp[7].split(":") #chr8:110351549:+
    chr1, pos1, strand1  = left_chr_bp_strand 
    chr2, pos2, strand2 = right_chr_bp_strand
    chr1 = chr1[3:]
    chr2 = chr2[3:]
    split_cnt = tmp[1]
    pair_cnt = tmp[2] 
    genes = tmp[0].split("--")
    gene1, gene2 = genes[0], genes[1]
    
#    to_write = [chr1, pos1, strand1, chr2, pos2, strand2, "RNA", args.sample_name, "NA", "NA", tool, split_cnt, pair_cnt, gene1, "NA", gene2, "NA"].join("\t")
    to_write = "\t".join([chr1, pos1, strand1, chr2, pos2, strand2, "RNA", "smc_rna_sim45", "NA", "NA", tool, split_cnt, pair_cnt, gene1, "NA", gene2, "NA"]) + "\n"
    cff_file.write(to_write) 
cff_file.close()
#    print >> out_file, "\t".join(map(str, [fusion.chr1, fusion.pos1, fusion.strand1, fusion.chr2, fusion.pos2, fusion.strand2, "RNA", sample, sample_type, disease_name, args.tool, fusion.split_cnt, fusion.pair_cnt, fusion.gene1, fusion.gene_location1, fusion.gene2, fusion.gene_location2]))

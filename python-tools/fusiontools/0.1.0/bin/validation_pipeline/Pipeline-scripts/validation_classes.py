#! /usr/bin/env python

class ValidatedFusion():
    """
    Represents one line of a validated fusion file. Takes in a line in the following format:
    gene1   gene2   strand1   strand2   chr1    chr2    start1  end1    start2  end2
    """
    def __init__(self, validated_fusion_line):
        self.__load_validated_fusion(validated_fusion_line)

    def __load_validated_fusion(self, validated_fusion_line):
        line = validated_fusion_line.split()
        #print(validated_fusion_line)
        self.gene1 = line[0]
        self.gene2 =line[1]
        self.gene1_candidates = ""
        self.gene2_candidates = ""
        self.gene1_strands = ""
        self.gene2_strands = ""
        #self.gene1_strands = line[2]
        #self.gene2_strands = line[3]
        #self.ensg1 = line[2]
        #self.ensg2 = line[3]
        #self.stands1 = line[2] 
        #self.strands2 = line[3]
        self.chr1 = line[4]
        self.chr2 = line[5]
        self.start1 = int(line[6])
        self.end1 = int(line[7])
        self.start2 = int(line[8])
        self.end2 = int(line[9])
        self.line = validated_fusion_line
        # attributes to collect all possible gene mappings from .bed file
        self.left = []
        self.right = []
        self.sample_name = line[10]
        if len(line) > 11:
            self.gene1_candidates = line[11]
            self.gene2_candidates = line[12]

class FusionInspectorFusion():
    """
    Represents one line of a FusionInspector output file finspector.fusion_predictions.final.abridged.FFPM.
    #FusionName     JunctionReadCount       SpanningFragCount       LeftGene        LeftLocalBreakpoint     LeftBreakpoint  RightGene       RightLocalBreakpoint    RightBreakpoint SpliceType      LargeAnchorSupport      NumCounterFusionLeft    NumCounterFusionRight   FAR_left        FAR_right       LeftBreakDinuc  LeftBreakEntropy        RightBreakDinuc RightBreakEntropy       FFPM
    """
    def __init__(self, FI_output_line):
        self.__load_FI_fusion(FI_output_line)

    def __load_FI_fusion(self, FI_output_line):
        self.line = FI_output_line
        line = FI_output_line.split()
        self.gene_pair = line[0]
        self.gene1 = line[0].split("--")[0]
        self.gene2 = line[0].split("--")[1]
        # chr17:61906923:+
        self.chr1 = line[5].split(":")[0]
        self.chr2 = line[8].split(":")[0]
        self.bp1 = int(line[5].split(":")[1])
        self.bp2 = int(line[8].split(":")[1])
        self.strand1 = line[5].split(":")[2]
        self.strand2 = line[8].split(":")[2]

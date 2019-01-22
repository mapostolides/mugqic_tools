#!/bin/bash

#    A script to convert cff input file to .cluster file format. Skips reannotation/filtering steps normally associated with .cluster files. Allows for direct examination of fusion caller results by using generated .cluster file as input to validate_fusion_stats.py
#    CFF column headers:
#    chr1    pos1    strand1 chr2    pos2    strand2 library sample_name     sample_type     disease tool    split_cnt       span_cnt        t_gene1 t_area1 t_gene2 t_area2
#    17      48539663        +       8       42620347        +       RNA     smc_rna_sim52   Tumor   VALIDATION      defuse  23      3       ACSF2   coding  CHRNA6  coding
#
#    CLUSTER column headers:
#    cluster_type    gene1   gene2  max_split_cnt max_span_cnt sample_type disease tools inferred_fusion_type gene1_on_bnd gene1_close_to_bnd  gene2_on_bnd gene2_close_to_bnd dna_supp  samples chr1 breakpoint_1 chr2 breakpoint_2
#    Gene_Cluster    LINC00359       LINC00359       39      17      Tumor   VALIDATION      defuse  SameGene        True    True    True    True    -1      smc_rna_sim52   chr13   97608116        chr13   97594053


awk '{print "NULL\t"$14"\t"$16"\t"$12"\t"$13"\t"$9"\t"$10"\t"$11"\t""NULL\t""NULL\t""NULL\t""NULL\t""NULL\t""-1\t"$8"\tchr"$1"\t"$2"\tchr"$4"\t"$5}' $1

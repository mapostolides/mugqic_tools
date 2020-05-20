import pandas as pd
import  pygeneann_reads_capture_DEV as pygeneann
import pybedtools.bedtool as bedtools
cff="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/outdir.TEST_FUSIONS.May-19-2020-WITH_RENAME/BCL3--CTB-171A8.1.cff.renamed.reann"
cff="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/testing_pipeline/outdir.TEST_FUSIONS.May-19-2020-WITH_RENAME/BCL3--CTB-171A8.1.cff.renamed.reann.MOD"
lines=[line for line in open(cff, "r")]
fusion=pygeneann.CffFusion(lines[0])
#header=['chr1', 'pos1', 'strand1', 'chr2', 'pos2', 'strand2', 'library', 'sample_name', 'sample_type', 'disease', 'tool', 'split_cnt', 'span_cnt', 't_gene1', 't_area1', 't_gene2', 't_area2', 'category', 'reann_gene1', 'reann_type1', 'reann_gene2', 'reann_type2', 'gene1_on_bdry', 'gene1_close_to_bndry', 'gene2_on_bdry', 'gene2_close_to_bndry', 'score', 'coding_id_distance', 'gene_interval_distance', 'dna_support', 'fusion_id', 'seq1', 'seq2', 'is_inframe', 'splice_site1', 'splice_site2', 'captured_reads']
header=fusion.zone1_attrs + fusion.zone2_attrs + fusion.zone3_attrs + fusion.zone4_attrs
#print(header)
#df=pd.read_csv(cff, sep='\t', keep_default_na=False, names=header)
df=pd.read_csv(cff, sep='\t', keep_default_na=False, index_col=False, names=header)
#print(df.to_string())
#create BedTools object with appropriate column names
df_bed=df[['chr1','pos1','pos1','chr2','pos2','pos2', 'fusion_id']]
df_bed.columns=['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id']
df_bed.loc[:,['pos1_2','pos2_2']] +=1
#print(df_bed.to_string())
df_bed=bedtools.BedTool.from_dataframe(df_bed)
#print(df_bed)
#print(bedtools.BedTool.from_dataframe(df_bed))
#string=fusion.tostring()

header = ['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id', 'chr1.1','pos1.1','pos1_2.1','chr2.1','pos2.1','pos2_2.1', 'fusion_id.1']
df_intersect=df_bed.pair_to_pair(df_bed, slop=100, rdn=True)
#print(df_bed.pair_to_pair(df_bed, slop=100, rdn=True))
df=df_intersect.to_dataframe(header=None)
df=df.iloc[:,0:14]
df.columns = header 
df_grouped=df.groupby('fusion_id')['fusion_id.1'].apply(lambda x: ','.join(x))
print(df_grouped[['fusion_id']]
print(df_grouped)
#print(df.to_string())


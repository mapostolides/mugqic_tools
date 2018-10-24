# load file, with LHS and RHS gene pairs
fusion_file <- read.table(file='integrate-validated-fusions.txt', sep = '\t', header= T, check.names=F)

left_genes <- as.character(fusion_file[[1]])
right_genes <- as.character(fusion_file[[2]])

library("biomaRt")
# use the hg19/grch37 database
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")

# get genomic start and end of a gene by using its hgnc_symbol
attributes <- c('chromosome_name', 'start_position', 'end_position', 'ensembl_gene_id', 'hgnc_symbol')
left_data<-getBM(attributes=attributes, filters="hgnc_symbol", values=left_genes, mart=ensembl)
right_data<-getBM(attributes=attributes, filters="hgnc_symbol", values=right_genes, mart=ensembl)

# get genes not found in database for both left and right data
left_diff <- setdiff(left_genes, left_data$hgnc_symbol)
right_diff <- setdiff(right_genes, right_data$hgnc_symbol)

left_diff_data <- getBM(attributes=attributes, filters="ensembl_gene_id", values=left_diff$ensID, mart=ensembl)
# order properly
left_diff_data <- left_diff_data[match(left_diff$ensID, left_diff_data$ensembl_gene_id),]

#generate data frame for final file
fusion_df <- data.frame(gene1=left_data$hgnc_symbol, gene2=right_data$hgnc_symbol, ensg1=left_data$ensembl_gene_id, ensg2=right_data$ensembl_gene_id, chr1=left_data$chromosome_name, chr2=right_data$chromosome_name, start1=left_data$start_position, end1=left_data$end_position, start2=right_data$start_position, end2=right_data$end_position)

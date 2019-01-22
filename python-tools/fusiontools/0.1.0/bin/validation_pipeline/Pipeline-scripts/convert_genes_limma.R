#! /Library/Frameworks/R.framework/Resources/bin/Rscript

library("limma")

# instantiate command line arguments
args <- commandArgs(TRUE)
file <- as.character(args[1])


alias2SymbolTable_null <-function(name){
	symbol<-alias2SymbolTable(name)
	if(length(symbol)==0){
		symbol<- name
	}
	return(symbol)
}

cff_file <- read.table(file=file, sep = '\t', header= F, check.names=F, stringsAsFactors=F)

left_genes <- unlist(cff_file[[14]])
right_genes <- unlist(cff_file[[16]])

official_names_left_genes <- alias2SymbolTable(left_genes)
official_names_right_genes <- alias2SymbolTable(right_genes)

# printing to stdout, gets read in by calling python script  
print( paste(official_names_left_genes, collapse=" " ) )
print( paste(official_names_right_genes, collapse=" " ) )

# match original list to official names list to account for NA values returned by alias2SymbolTable 


# SCRAP


#alias2SymbolTablev <- Vectorize(alias2SymbolTable)
#names=alias2SymbolTablev(alias_vector)
#print(names)
#alias_vs_name<-data.frame(alias=alias_vector, names)
#unlist(lapply(cff_file[14], alias2SymbolTable_null))
#print(cff_file[14])
#print(lapply(cff_file[14], alias2SymbolTable))


#for (name in cff_file[14]){
#	name_orig <- name
#	print(paste("before", name, sep=" "))
#	symbol<-alias2SymbolTable(name)
#	print(paste("after", symbol, sep=" "))

	#print(paste(symbol, name_orig, sep=" "))
#}

#for (name in cff_file[14])  print(alias2SymbolTable_null(name))
#alias2SymbolTable_null("EMR3")
#mclapply(names, wrapper, cores=threads)



#cff_file[16]

#left_genes
#right_genes


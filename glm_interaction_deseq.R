#given a set of genes with replicates for both condition and genotype
#calls differential expression, conditioning purely on condition and genotype
#calls differential expression based on condition x genotype interaction, via likelihood ratio test
#calls fold change for each of the above comparisons
#generates volcano plot of each of the above comparisons

#primary dependency is bioconductor, r package
#additional requirement, for volcano plot, is maDB

library ( "DESeq" ) 

#specify file with one gene per row and counts per condition in each column
countfile = "genomecounts.txt" 

#read in the table specifying the counts for each gene in each condition
atxaCountTable <- read.table ( countfile, header=TRUE,row.names=1)

#specify condition and genotype assignments for each sample
atxaDesign <- data.frame( row.names = colnames( atxaCountTable ), condition = c( "CO2", "CO2", "air", "air", "air", "air", "CO2", "CO2" ), genotype = c( "Ames35", "Ames35", "dAtxA", "dAtxA", "Ames35", "Ames35", "dAtxA", "dAtxA" ) )

#specify design matrix of experiment
cdsFull <- newCountDataSet ( atxaCountTable, atxaDesign )

cdsFull <- estimateSizeFactors ( cdsFull )

cdsFull <- estimateDispersions ( cdsFull ) 

#function to plot (log) dispersion estimates versus (log) read count means
plotDispEsts <- function (cdsFull ) { 
plot( rowMeans(counts (cdsFull, normalized=TRUE ) ), fitInfo(cdsFull)$perGeneDispEsts, pch = '.', log="xy" ) 
#creates best fit line for visualization, though these are not employed in inference if sharingMode = "maximum"
xg <- 10^seq( -.5, 5, length.out=300 ) 
lines( xg, fitInfo(cdsFull)$dispFun( xg ), col="red" )
}

#un-comment to call the plot
#plotDispEsts ( cdsFull )

fit0a <- fitNbinomGLMs( cdsFull, count ~ genotype ) 
fit0b <- fitNbinomGLMs( cdsFull, count ~ condition ) 
fit1 <- fitNbinomGLMs( cdsFull, count ~ genotype + condition )
fit2 <- fitNbinomGLMs( cdsFull, count ~ genotype + condition + genotype:condition )

#compare fit1 to fit0b = effect of condition
pvalsGLM_condition <- nbinomGLMTest( fit1, fit0a )
#control false discovery rate with benjamini-hochberg
padjGLM_condition <- p.adjust( pvalsGLM_condition, method="BH") 
#compare fit1 to fit0a = effect of genotype
pvalsGLM_genotype <- nbinomGLMTest( fit1, fit0b )
#control false discovery rate with benjamini-hochberg
padjGLM_genotype <- p.adjust( pvalsGLM_genotype, method="BH") 
#compare #fit2 to fit1 = effect of interaction
pvalsGLM_interaction <- nbinomGLMTest( fit2, fit1 )
#control false discovery rate with benjamini-hochberg
padjGLM_interaction  <- p.adjust( pvalsGLM_interaction, method="BH") 

#generates volcano plot for condition 
drawVolcanoPlot(fit2$conditionair,padjGLM_condition)
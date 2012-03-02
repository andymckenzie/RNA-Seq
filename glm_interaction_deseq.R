#given a set of genes with replicates for both condition and genotype
#calls differential expression, conditioning purely on condition and genotype
#calls differential expression based on condition x genotype interaction, via likelihood ratio test
#calls fold change for each of the above comparisons
#generates volcano plot of each of the above comparisons

#primary dependency is bioconductor, r package

library ( "DESeq" ) 
library ( ggplot2 ) 

#free parameter
expression_cutoff = 0.3

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

#filter out genes with low expression from both the data frame and the table
rs <- rowSums ( counts ( cdsFull )) 
use <- (rs > quantile(rs, expression_cutoff)) 

#filter genes to remove those less than at the 0.3 percentile of total row sum count
#to a good approx, this measure is uncorrelated with the test statistics under the null hypothesis
cdsFilt <- cdsFull [ use, ]
rowsum_count = rowSums (atxaCountTable) 
use_count = (rowsum_count > quantile(rowsum_count, expression_cutoff))
atxaCountTable_Filt = atxaCountTable [ use_count, ]  

fit0a <- fitNbinomGLMs( cdsFilt, count ~ genotype ) 
fit0b <- fitNbinomGLMs( cdsFilt, count ~ condition ) 
fit1 <- fitNbinomGLMs( cdsFilt, count ~ genotype + condition )
fit2 <- fitNbinomGLMs( cdsFilt, count ~ genotype + condition + genotype:condition )

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

#calculate fold change in Ames35_CO2 vs dAtxa_air
Ames35_CO2_mean = apply(cbind(atxaCountTable_Filt$Ames35_CO2_1, atxaCountTable_Filt$Ames35_CO2_2), 1, mean) 
dAtxA_air_mean = apply(cbind(atxaCountTable_Filt$Ames35delAtxA_Air_1, atxaCountTable_Filt$Ames35delAtxA_Air_2), 1, mean)
fold_change = Ames35_CO2_mean / dAtxA_air_mean

#highlight genes that are in the "extended pathogenicity island" 
atxaCountTable_Filt$threshold = as.factor(row.names(atxaCountTable_Filt) > "RBAH05429" & row.names(atxaCountTable_Filt) < "RBAH05492")
 
#construct the "volcano" plot 
table = cbind(padjGLM_interaction, fold_change) 
g = ggplot(data=atxaCountTable_Filt, aes(x=log2(fold_change), y=-log10(padjGLM_interaction), colour=threshold)) +
  geom_point(alpha=0.4, size=2.5) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) + 
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw() + opts(legend.position="none") 
g
library(ggplot2)

gf = read.table("bayes_dstn_gene_of_interest.txt", header = F)

#str(gf), yields 'data,frame' : 400 obs of 4 variables, $distance : int 1 2 3 4 5 ... $prior num 0.000843 , etc for $likelihood, $normpost

UTR_Length = gf$V1

expprior = (gf$V2)
Prior = expprior / sum(expprior)

explikelihood = (gf$V3)
Likelihood = explikelihood / sum(explikelihood)

expposterior = (gf$V4)
Posterior = expposterior / sum(expposterior)

gh = cbind(UTR_Length, Prior, Likelihood, Posterior) 

gg = data.frame(gh)

dfm = melt (gg, id = "UTR_Length", measure = c("Prior", "Likelihood", "Posterior"))

a = qplot(UTR_Length, value, data = dfm, geom = "line", colour = variable) + scale_x_continuous('Upstream Untranslated Region Length',breaks = c(0,100,200,300,400)) + scale_y_continuous('Probability Of Being TSS')
library(gld)

utr = read.table("utr_length.txt")

fit.result <- starship(utr$V2)

lambda1 <- fit.result$optim.results$par[1]
lambda2 <- fit.result$optim.results$par[2]
lambda3 <- fit.result$optim.results$par[3]
lambda4 <- fit.result$optim.results$par[4]

norm_utr = utr$V2 / sum(utr$V2)

plot(x = utr$V1, y = norm_utr, ylab = "Probability Of Being TSS", xlab = "5' Untranslated Region Length", xlim =c(0,400))
#add a scale parameter theta = 1.25 to adjust for better fit, visually
curve(1.25*dgl(x/5,lambda1,lambda2,lambda3,lambda4),add=TRUE)

utrseq = 1:400

prior = (1.25*dgl((utrseq/5),lambda1,lambda2,lambda3,lambda4))

logprior = log(prior)

write(logprior, file = "bayes_gld_tss_log_prior.txt", ncolumns = 1)
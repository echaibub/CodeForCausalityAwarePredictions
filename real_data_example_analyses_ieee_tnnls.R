####################################################
## run real data experiments
####################################################

## The data used for the the real data illustrations is available at: 
## https://www.synapse.org/#!Synapse:syn35180598/files/ 
##  
## This dataset represents a subset of the entire mPower data, and was 
## collected from human subjects who consented to share the data with 
## qualified researchers. The mPower data is publicly available to anyone 
## who agrees to specific conditions of use and goes  through a 
## qualification process described in detail at: 
## https://www.synapse.org/#!Synapse:syn4993293/wiki/247860. 


source("utility_functions_ieee_tnnls.R")

require(synapser)
synLogin()

load(synGet("syn35180808")$path)

dat$dage <- GetDiscretizedAge(dat, n.levels = 3, breaks = c(17, 44, 65, 99), 
                              level.names = c("YoungAge", "MiddleAge", "SeniorAge"))


nRuns <- 100

set.seed(123456789)
my.seeds <- sample(c(10000:100000), nRuns, replace = FALSE)

feature.names <- names(dat)[4:11]

output.path <- ""

file.name <- paste(output.path, "output_tnnls_none.RData", sep = "")
out <- RunAnalyses(dat, 
                   label.name = "PD", 
                   feature.names = feature.names, 
                   confounder.name = "dage",
                   neg.class.name = "FALSE", 
                   pos.class.name = "TRUE",
                   run.seeds = my.seeds,
                   adjustment.method = "none",
                   prop = 0.8)
save(out, file = file.name, compress = TRUE)


file.name <- paste(output.path, "output_tnnls_causality_aware_lm.RData", sep = "")
out <- RunAnalyses(dat, 
                   label.name = "PD", 
                   feature.names = feature.names, 
                   confounder.name = "dage",
                   neg.class.name = "FALSE", 
                   pos.class.name = "TRUE",
                   run.seeds = my.seeds,
                   adjustment.method = "causality.aware.lm",
                   prop = 0.8)
save(out, file = file.name, compress = TRUE)


## linear model baseline 
file.name <- paste(output.path, "output_tnnls_causality_aware_lm_train_only.RData", sep = "")
out <- RunAnalyses(dat, 
                   label.name = "PD", 
                   feature.names = feature.names, 
                   confounder.name = "dage",
                   neg.class.name = "FALSE", 
                   pos.class.name = "TRUE",
                   run.seeds = my.seeds,
                   adjustment.method = "causality.aware.lm.train.only",
                   prop = 0.8)
save(out, file = file.name, compress = TRUE)


file.name <- paste(output.path, "output_tnnls_causality_aware_gam.RData", sep = "")
out <- RunAnalyses(dat, 
                   label.name = "PD", 
                   feature.names = feature.names, 
                   confounder.name = "age",
                   neg.class.name = "FALSE", 
                   pos.class.name = "TRUE",
                   run.seeds = my.seeds,
                   adjustment.method = "causality.aware.gam",
                   prop = 0.8)
save(out, file = file.name, compress = TRUE)


## additive model baseline 
file.name <- paste(output.path, "output_tnnls_causality_aware_gam_train_only.RData", sep = "")
out <- RunAnalyses(dat, 
                   label.name = "PD", 
                   feature.names = feature.names, 
                   confounder.name = "age",
                   neg.class.name = "FALSE", 
                   pos.class.name = "TRUE",
                   run.seeds = my.seeds,
                   adjustment.method = "causality.aware.gam.train.only",
                   prop = 0.8)
save(out, file = file.name, compress = TRUE)


file.name <- paste(output.path, "output_tnnls_matching.RData", sep = "")
out <- RunAnalyses(dat, 
                   label.name = "PD", 
                   feature.names = feature.names, 
                   confounder.name = "dage",
                   neg.class.name = "FALSE", 
                   pos.class.name = "TRUE",
                   run.seeds = my.seeds,
                   adjustment.method = "matching",
                   prop = 0.8)
save(out, file = file.name, compress = TRUE)


file.name <- paste(output.path, "output_tnnls_approx_IPW_ps.RData", sep = "")
out <- RunAnalyses(dat, 
                   label.name = "PD", 
                   feature.names = feature.names, 
                   confounder.name = "dage",
                   neg.class.name = "FALSE", 
                   pos.class.name = "TRUE",
                   run.seeds = my.seeds,
                   adjustment.method = "approximate.IPW",
                   prop = 0.8)
save(out, file = file.name, compress = TRUE)


#############################################
## generate figures
#############################################

file.names <- paste(output.path, c("output_tnnls_none.RData",
                                   "output_tnnls_causality_aware_lm_train_only.RData",
                                   "output_tnnls_causality_aware_gam_train_only.RData",
                                   "output_tnnls_matching.RData",
                                   "output_tnnls_approx_IPW_ps.RData",
                                   "output_tnnls_causality_aware_lm.RData",
                                   "output_tnnls_causality_aware_gam.RData"), sep = "")


method.names <- c("no adjustment",
                  "lin. model baseline",
                  "add. model baseline",
                  "matching",
                  "approximate IPW",
                  "lin. m. causality-aware",
                  "add. m. causality-aware")

n.runs <- 100

aucs.1 <- CatenateAucs1(file.names, n_runs = n.runs, method.names)
aucs.2 <- CatenateAucs2(file.names, n_runs = n.runs, method.names)

n.adj <- ncol(aucs.1)

fig.path <- ""


jit <- 0.1
my.cex <- 0.5
#pdf(paste(fig.path, "fig17.pdf", sep = ""), width = 3.0, height = 4.6)
par(mfrow = c(1, 1), mar = c(10, 3.5, 1, 0.5), mgp = c(2.5, 0.75, 0))
boxplot(aucs.1, las = 2, border = rgb(0, 0.5, 0.5, 1), col = "white", ylim = c(0.65, 0.9),
        main = "", ylab = "AUROC", xaxt = "n", at = seq(n.adj)-jit, outline = FALSE)
boxplot(aucs.2, las = 2, border = rgb(0.5, 0.5, 0, 1), col = rgb(0.5, 0.5, 0, 0.25), at = seq(n.adj)+jit, 
        add = TRUE, xaxt = "n", outline = FALSE)
axis(side = 1, at = seq(n.adj), labels = colnames(aucs.1), las = 2)
legend("topright", legend = c("no shift test set", "shifted test set"), 
       text.col = c(rgb(0, 0.5, 0.5, 1), rgb(0.5, 0.5, 0, 1)), bty = "n", cex = 1.0,
       fill = c("white", rgb(0.5, 0.5, 0, 0.25)),
       border = c(rgb(0, 0.5, 0.5, 1), rgb(0.5, 0.5, 0, 1)))
#dev.off()



dat2 <- dat
dat2$PD <- as.character(dat2$PD)
dat2$PD[dat$PD == "TRUE"] <- "PD"
dat2$PD[dat$PD == "FALSE"] <- "non-PD"
split.dat <- TrainTest1Test2DataSplit(dat = dat2, 
                                      label.name = "PD",
                                      confounder.name = "dage",
                                      neg.class.name = "non-PD",
                                      pos.class.name = "PD",
                                      prop = 0.8)

idx.train <- split.dat$idx.train
idx.test.1 <- split.dat$idx.test.1
idx.test.2 <- split.dat$idx.test.2


#pdf(paste(fig.path, "fig16.pdf", sep = ""), width = 6, height = 2.5)
par(mfrow = c(1, 3), mar = c(4, 3, 2, 0.5))
mosaicplot(PD ~ dage, data = dat2[idx.train,], las = 2, color = TRUE,
           main = "training set", ylab = "discretized age", xlab = "PD status")
mtext(side = 3, "(a)", at = 0)
mosaicplot(PD ~ dage, data = dat2[idx.test.1,], las = 2, color = TRUE, 
           main = "no shift test set", ylab = "discretized age", xlab = "PD status")
mtext(side = 3, "(b)", at = 0)
mosaicplot(PD ~ dage, data = dat2[idx.test.2,], las = 2, color = TRUE, 
           main = "shifted test set", ylab = "discretized age", xlab = "PD status")
mtext(side = 3, "(c)", at = 0)
#dev.off()


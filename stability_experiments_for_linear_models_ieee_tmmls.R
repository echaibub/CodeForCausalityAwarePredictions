source("utility_functions_ieee_tnnls.R")


###################################################
## experiments under the method's assumptions
###################################################

n <- 1e+3
n.sim <- 1e+3 
n.test <- 9
n.feat <- 5

gamma.seq <- c(seq(0, 1, length.out = 50), seq(2, 100, by = 2))

## fixed Var(Y_ts) experiments
ex1 <- RunShiftExperiments1(n,
                            n.sim,
                            n.feat,
                            response.name = "Y", 
                            confounder.name = "C",
                            var.C.vals = seq(1, 3, length.out = n.test), 
                            var.Y.vals = rep(1, n.test), 
                            cov.YC.vals = seq(0.8, -0.8, length.out = n.test),
                            beta.XC.range = c(0, 1),
                            beta.XY.range = c(0, 1),
                            rho.X.range = c(-0.8, 0.8),
                            my.seed = 123456789,
                            gamma.seq = gamma.seq)


ex2 <- RunShiftExperiments1(n,
                            n.sim,
                            n.feat,
                            response.name = "Y", 
                            confounder.name = "C",
                            var.C.vals = seq(1, 3, length.out = n.test), 
                            var.Y.vals = rep(1, n.test), 
                            cov.YC.vals = seq(0.8, -0.8, length.out = n.test),
                            beta.XC.range = c(-1, 0),
                            beta.XY.range = c(0, 1),
                            rho.X.range = c(-0.8, 0.8),
                            my.seed = 123456789,
                            gamma.seq = gamma.seq)


ex3 <- RunShiftExperiments1(n,
                            n.sim,
                            n.feat,
                            response.name = "Y", 
                            confounder.name = "C",
                            var.C.vals = rep(1, n.test), 
                            var.Y.vals = seq(1, 3, length.out = n.test), 
                            cov.YC.vals = seq(0.8, -0.8, length.out = n.test),
                            beta.XC.range = c(0, 1),
                            beta.XY.range = c(0, 1),
                            rho.X.range = c(-0.8, 0.8),
                            my.seed = 123456789,
                            gamma.seq = gamma.seq)


## increasing Var(Y_ts) experiments
ex4 <- RunShiftExperiments1(n,
                            n.sim,
                            n.feat,
                            response.name = "Y", 
                            confounder.name = "C",
                            var.C.vals = rep(1, n.test), 
                            var.Y.vals = seq(1, 3, length.out = n.test), 
                            cov.YC.vals = seq(0.8, -0.8, length.out = n.test),
                            beta.XC.range = c(-1, 0),
                            beta.XY.range = c(0, 1),
                            rho.X.range = c(-0.8, 0.8),
                            my.seed = 123456789,
                            gamma.seq = gamma.seq)


fig.path <- ""

fnm <- paste0(fig.path, "fig7.pdf")
#pdf(file = fnm, width = 6, height = 3.75)
MyPlot1(e = ex1, my.ylim.1 = c(0, 2.2), my.ylim.2 = c(0, 1))
#dev.off()

fnm <- paste0(fig.path, "fig8.pdf")
#pdf(file = fnm, width = 6, height = 3.75)
MyPlot1(e = ex2, my.ylim.1 = c(0, 2), my.ylim.2 = c(0, 0.5))
#dev.off()

fnm <- paste0(fig.path, "fig9.pdf")
#pdf(file = fnm, width = 6, height = 3.75)
MyPlot1(e = ex3, my.ylim.1 = c(0, 2), my.ylim.2 = c(0, 1))
#dev.off()

fnm <- paste0(fig.path, "fig10.pdf")
#pdf(file = fnm, width = 6, height = 3.75)
MyPlot1(e = ex4, my.ylim.1 = c(0, 2.5), my.ylim.2 = c(0, 1.1))
#dev.off()



###################################################
## failure mode experiments
###################################################

n <- 1e+3
n.sim <- 1e+3 
n.test <- 9

gamma.seq <- c(seq(0, 1, length.out = 50), seq(2, 100, by = 2))
beta.XC.test.values <- seq(1, 3, length.out = n.test)


## fixed Var(Y_ts) experiments
fm1 <- RunShiftExperimentsFailureMode(n,
                                      n.sim,
                                      response.name = "Y", 
                                      confounder.name = "C",
                                      beta.YC.range = c(0, 1),
                                      beta.XC.range = c(0, 1),
                                      beta.XY.range = c(0, 1),
                                      my.seed = 123456789,
                                      gamma.seq = gamma.seq,
                                      beta.XC.test.values)



fm2 <- RunShiftExperimentsFailureMode(n,
                                      n.sim,
                                      response.name = "Y", 
                                      confounder.name = "C",
                                      beta.YC.range = c(-1, 0),
                                      beta.XC.range = c(0, 1),
                                      beta.XY.range = c(0, 1),
                                      my.seed = 123456789,
                                      gamma.seq = gamma.seq,
                                      beta.XC.test.values)



fig.path <- ""

fnm <- paste0(fig.path, "fig12.pdf")
#pdf(file = fnm, width = 6, height = 3.75)
MyPlot1(e = fm1, my.ylim.1 = c(0.5, 2.5), my.ylim.2 = c(0, 1.2))
#dev.off()

fnm <- paste0(fig.path, "fig13.pdf")
#pdf(file = fnm, width = 6, height = 3.75)
MyPlot1(e = fm2, my.ylim.1 = c(0.5, 5.5), my.ylim.2 = c(0, 2))
#dev.off()



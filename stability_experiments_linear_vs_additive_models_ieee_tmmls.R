#####################################################
## run experiments
#####################################################

source("utility_functions_ieee_tnnls.R")

n <- 1e+3
n.runs <- 1e+3
n.feat <- 5

ae1 <- RunShiftExperiments2(n,
                            n.runs,
                            n.feat,
                            rho.YC.train = 0.8,
                            rho.YC.test = seq(0.8, -0.8, by = -0.2),
                            beta.XC.range = c(1, 3),
                            beta.XY.range = c(1, 3),
                            generate.non.linear.data = TRUE,
                            my.seed = 987654321)


ae2 <- RunShiftExperiments2(n,
                            n.runs,
                            n.feat,
                            rho.YC.train = 0.8,
                            rho.YC.test = seq(0.8, -0.8, by = -0.2),
                            beta.XC.range = c(-3, -1),
                            beta.XY.range = c(1, 3),
                            generate.non.linear.data = TRUE,
                            my.seed = 123456789)


#####################################
## generate figures
#####################################

fig.path <- ""

fnm <- paste0(fig.path, "fig14.pdf")
#pdf(file = fnm, width = 7, height = 3.0)
MyPlot2(MSEs.lm = ae1$MSEs.lm,
      MSEs.lm.ca = ae1$MSEs.lm.ca,
      MSEs.gam = ae1$MSEs.gam,
      MSEs.gam.ca = ae1$MSEs.gam.ca,
      ylim1 = c(0.21, 1.3),
      ylim2 = c(0, 0.4))
#dev.off()


fnm <- paste0(fig.path, "fig15.pdf")
#pdf(file = fnm, width = 7, height = 3.0)
MyPlot2(MSEs.lm = ae2$MSEs.lm,
      MSEs.lm.ca = ae2$MSEs.lm.ca,
      MSEs.gam = ae2$MSEs.gam,
      MSEs.gam.ca = ae2$MSEs.gam.ca,
      ylim1 = c(0.21, 1.3),
      ylim2 = c(0, 0.4))
#dev.off()


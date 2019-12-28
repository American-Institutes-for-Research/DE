require(testthat)
set.seed(142857)
n <- 2000
theta <- rnorm(n)
x1 <- runif(n)
theta <- theta + x1 * 0.2
mean(theta)
gen <- function(paramTab, theta, key) {
  params <- DE:::csv2paramList(pModel='GPCM', paramTab, key)
  # get parameter values from items
  a <- params$'3pl'$a
  b <- params$'3pl'$b
  c <- params$'3pl'$c
  t(sapply(theta, function(thetai) {
    rbinom(length(a),1, c + (1 - c) / (1 + exp(-1.7 * a * (thetai - b))))
  }))
}
# paramTab
paramTab <- structure(list(ItemID = structure(1:13,
                                              .Label = c("m017401", "m017701", "m017901", "m018201", "m018401", "m018501", "m018601", "m020001", "m020501", "m046301", "m046501", "m051501", "n202831"), class = "factor"),
                      P0 = c(0.25, 1, 1.15, 0.52, 1.11, 1.64, 0.78, 0.72, 0.72, 0.89, 0.92, 1.2, 0.75),
                      P1 = c(-5.16, -1.01, -0.93, -1.21, -1.03, 0.34, 0.9, -0.49, -0.62, -1.07, -0.23, 1.22, -2.58),
                      P2 = c(0.19, 0.16, 0.15, 0.03, 0.24, 0.26, 0.12, 0, 0, 0.28, 0.33, 0.2, 0.25),
                      P3 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                      P4 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                      P5 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                      P6 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                      P7 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                      P8 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                      P9 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), 
                      P10 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                      ScorePoints = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
                      MODEL = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "3pl", class = "factor")), class = "data.frame", row.names = c(NA, -13L))
mat <- gen(paramTab, theta, key=paramTab$ItemID)
colnames(mat) <- c("m017401", "m017701", "m017901", "m018201", "m018401", "m018501", 
"m018601", "m020001", "m020501", "m046301", "m046501", "m051501", 
"n202831")
mat <- data.frame(mat)
rownames(mat) <- paste0("pseudo-student",1:nrow(mat))
mat$origwt <- 1
nperstratum <- 10
nstrata <- length(theta)/nperstratum
mat$repgrp1 <- rep(1:nstrata, each=nperstratum)
mat$jkunit <- rep(rep(1:2, each=nperstratum/2), nstrata)
stuItems <- mat[,1:13]
stuItems$oppID <- factor(rownames(mat), levels=rownames(mat))
stuItems <- reshape(data=stuItems, varying=c(paramTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=paramTab$ItemID, timevar="key")
rownames(stuItems) <- NULL
stuDat <- mat[, c('origwt', 'repgrp1', 'jkunit')]
stuDat$oppID <- rownames(stuDat)
############### test functions ###############
attr(paramTab, "location") <- 277.1563
attr(paramTab, "scale") <- 37.7297
mat$x1 <- stuDat$x1 <- x1
stuDat$origwt <- mat$origwt <- runif(nrow(stuDat)) * 4 * abs(stuDat$x1 + 3)

# tests:
context("mml works") #When this fails all regression tests are invalid.
mml1 <- mml(stuItems=split(stuItems, stuItems$oppID), stuDat=stuDat, paramTab=paramTab, Q=34, idVar="oppID")
expect_is(mml1, "mml.means")
sis <- split(stuItems, stuItems$oppID)
context("mml errors") #When this fails all regression tests are invalid.
expect_error(mml1 <- mml(stuItems=sis[1:10], stuDat=stuDat, paramTab=paramTab, regType = 'regression', Q=34, idVar="oppID"), "pseudo-student11")
expect_error(mml1 <- mml(stuItems=sis, stuDat=stuDat[-11,], paramTab=paramTab, regType = 'regression', Q=34, idVar="oppID"), "pseudo-student11")

context("mml summary and variance estimation") #When this fails all regression tests are invalid.
expect_equal(coef(mml1), c(`(Intercept)` = 0.116721240655054, "Population SD" = 0.972350274042811))
# gradientHessian is a quirk of AM, set it to TRUE to compare the apples to apples
mml1s <- summary(mml1, gradientHessian=TRUE)
expect_is(mml1s, "summary.mml.means")
# values from AM with convergence set to 1E-12
mml1s_coef_REF <- structure(c(281.560163073, 36.6864807567,
                              0.974354739545, 0.93732233951),
                            .Dim = c(2L,2L), .Dimnames = list(c("(Intercept)", "Population SD"),
                                                         c("Estimate", "StdErr")))
# compare means
expect_equal(mml1s$coefficients[,1], mml1s_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*20)
# compare var estimates
expect_equal(mml1s$coefficients[,2], mml1s_coef_REF[,2], tolerance=sqrt(.Machine$double.eps)*20)

###
mml1R_SE_REF <- c(`(Intercept)` =  0.97801507512794, `Population SD` = 0.916317157364324)
mml1Robust <- summary(mml1, varType="robust")
expect_equal(mml1Robust$coef[,"StdErr"], mml1R_SE_REF)

# REF values from AM
mml1T_SE_REF <- c(`(Intercept)` = 0.967524019094, `Population SD` = 0.873102164238)
mml1Taylor <- summary(mml1, varType="Taylor", stratavar="repgrp1", psuvar="jkunit")
expect_equal(mml1Taylor$coef[,"StdErr"], mml1T_SE_REF, tolerance=0.02)

mml1C_SE_REF <- c(`(Intercept)` = 0.958703029704291, `Population SD` = 0.881032943695999)
mml1Cluster <- summary(mml1, varType="cluster", clustervar="repgrp1")
expect_equal(mml1Cluster$coef[,"StdErr"], mml1C_SE_REF)

mml2 <- mml(~x1, stuItems=split(stuItems, stuItems$oppID), stuDat=stuDat, paramTab=paramTab, regType = 'regression', Q=34, idVar="oppID")
expect_is(mml2, "mml.means")
expect_equal(coef(mml2), c(`(Intercept)` = 0.00592827512675005, x1 = 0.224670345051881, "Population SD" = 0.970126544720125))
mml2s <- summary(mml2, gradientHessian=TRUE)
# results form AM
mml2s_coef_REF <- structure(c(277.379973218, 8.47675551939, 36.602582391, 
                              1.91227846027, 3.42771590088, 0.938701562756),
                            .Dim = c(3L, 2L), .Dimnames = list(c("(Intercept)", "x1", "Population SD"),
                                                               c("Estimate", "StdErr")))
expect_equal(mml2s$coefficients[,1], mml2s_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*10)
expect_equal(mml2s$coefficients[,2], mml2s_coef_REF[,2], tolerance=sqrt(.Machine$double.eps)*10)
mml2R_SE_REF <- c(`(Intercept)` = 3.81354746413122, x1 = 6.53149834137551, `Population SD` = 0.9125654117294)
mml2Robust <- summary(mml2, varType="robust")
expect_equal(mml2Robust$coef[,"StdErr"], mml2R_SE_REF)

# gradientHessian=TRUE is what AM uses
mml2Taylor <- summary(mml2, varType="Taylor", stratavar="repgrp1", psuvar="jkunit", gradientHessian=TRUE)
mml2T_SE_REF <- c(`(Intercept)` = 1.91342587986, x1 = 3.35043639323, `Population SD` = 0.872891128653)
expect_equal(mml2Taylor$coef[,"StdErr"], mml2T_SE_REF, tolerance=0.001)
mml2Cluster <- summary(mml2, varType="cluster", clustervar="repgrp1")
mml2C_SE_REF <- c(`(Intercept)` = 1.88505689357269, x1 = 3.33275071678814, `Population SD` = 0.87139149522642)
expect_equal(mml2Cluster$coef[,"StdErr"], mml2C_SE_REF)

context("weighted case")

mml2W <- mml(~x1,stuItems=split(stuItems, stuItems$oppID), stuDat=stuDat, paramTab=paramTab, Q=34, idVar="oppID", weightvar="origwt")
mml2s_coef_REF <- structure(c(276.78360779, 9.90868559817, 36.9630033481, 
                              2.1549072924, 3.85424777114, 1.06630035114),
                            .Dim = c(3L, 2L), .Dimnames = list(c("(Intercept)", "x1", "Population SD"),
                                                               c("Estimate", "StdErr")))
mml2WTaylor <- summary(mml2W, varType="Taylor", stratavar="repgrp1", psuvar="jkunit", gradientHessian=TRUE)
expect_equal(mml2WTaylor$coefficients[,1], mml2s_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*20)
# compare var estimates
expect_equal(mml2WTaylor$coefficients[,2], mml2s_coef_REF[,2], tolerance=2*(.Machine$double.eps)^0.25)


context("factor")
b <- table(stuItems$oppID, stuItems$score)
stuDat$multiLevel <- factor(sample(c(1:5,NA), nrow(stuDat), prob=c(0.5,0.3,0.10,0.03,0.02,0.05), replace=TRUE), 1:5, LETTERS[1:5])
stuDat$mlM <- mat$mlM <- 0 + stuDat$multiLevel %in% NA
stuDat$mlE <- mat$mlE <- 0 + (stuDat$multiLevel %in% "E" | stuDat$oppID %in% rownames(b[b[,2] == 13,])[1:80])
stuDat$mlA <- mat$mlA <- 0 + (stuDat$multiLevel %in% "A" & !stuDat$mlE)
stuDat$mlB <- mat$mlB <- 0 + (stuDat$multiLevel %in% "B" & !stuDat$mlE)
stuDat$mlC <- mat$mlC <- 0 + (stuDat$multiLevel %in% "C" & !stuDat$mlE)
stuDat$mlD <- mat$mlD <- 0 + (stuDat$multiLevel %in% "D" & !stuDat$mlE)
stuDat$mlM <- mat$mlM <- 0 + (stuDat$multiLevel %in% NA & !stuDat$mlE)
stuItemsFactor <- stuItems
naitems <- c("m017401", "m017701", "m017901", "m018201", "m018401", "m018501", "m018601", "m020001", "m020501", "m046301", "m046501")
stuItemsFactor$score[stuItems$oppID %in% rownames(b[b[,2] == 13,])[1:20] & stuItems$key %in% naitems] <- NA
matFactor <- mat
for(nai in 1:length(naitems)) {
  matFactor[rownames(matFactor) %in% rownames(b[b[,2] == 13,])[1:20],naitems[nai]] <- NA
}

mmlml <- mml(~mlA + mlB + mlC + mlD + mlE, stuItems=split(stuItemsFactor, stuItems$oppID), stuDat=stuDat, paramTab=paramTab, Q=34, idVar="oppID", weight="origwt")
mmlmlsTaylor <- summary(mmlml, varType="Taylor", stratavar="repgrp1", psuvar="jkunit", gradientHessian=TRUE)
mmlmls_coef_REF <- structure(c(284.637143, -5.9554397, -6.52445584, -4.22802448, -10.1810361,  61.2638569, 34.9988724,
                               5.73766411, 5.90805667,  6.12720888,  6.63789466,  8.14112037,  7.30977647, 1.04592348),
                             .Dim = c(7L, 2L), .Dimnames = list(c("(Intercept)", "mlA", "mlB", "mlC", "mlD", "mlE", "Population SD"),
                                                                c("Estimate", "StdErr")))
 
expect_equal(mmlmlsTaylor$coefficients[,1], mmlmls_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*200)
expect_equal(mmlmlsTaylor$coefficients[,2], mmlmls_coef_REF[,2], tolerance=2*(.Machine$double.eps)^0.25)

context("missing data elements")

mat1 <- mat
mat1$m017401[1:1000] <- NA
stuItems2 <- mat1[,1:13]
stuItems2$oppID <- factor(rownames(mat1), levels=rownames(mat1))
stuItems2 <- reshape(data=stuItems2, varying=c(paramTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=paramTab$ItemID, timevar="key")
#stuItems2 <- gather(stuItems2,key,score, m017401:n202831) 

mml1B <- mml(stuItems=split(stuItems2, stuItems2$oppID), stuDat=stuDat, paramTab=paramTab, Q=34, idVar="oppID")
mml1Bs <- summary(mml1B, gradientHessian=TRUE)

# AM results with NA
mml1Bs_coef_REF <- structure(c(281.55602373, 36.6399666578,
                               0.973407969332, 0.940574435195),
                             .Dim = c(2L,2L), .Dimnames = list(c("(Intercept)", "Population SD"),
                                                          c("Estimate", "StdErr")))

expect_equal(mml1Bs$coefficients[,2], mml1Bs_coef_REF[,2], tolerance=20*sqrt(.Machine$double.eps))

# AM Taylor results with NA
mml1BsT_coef_REF <- structure(c(281.55602373, 36.6399666578,
                                0.96547714892, 0.877505088603),
                              .Dim = c(2L,2L), .Dimnames = list(c("(Intercept)", "Population SD"),
                                                                c("Estimate", "StdErr")))

mml1BTaylor <- summary(mml1B, varType="Taylor", stratavar="repgrp1", psuvar="jkunit", gradientHessian=TRUE)
expect_equal(mml1BTaylor$coefficients[,1], mml1BsT_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*20)
expect_equal(mml1BTaylor$coefficients[,2], mml1BsT_coef_REF[,2], tolerance=10*(.Machine$double.eps)^0.25)

context("missing data coded 8")

mat8 <- mat
set.seed(142857) # make sure sample is the same
randomRow <- sample(1:2000,1000,replace=TRUE) 
randomCol <- sample(1:13,1000,replace=TRUE)
# randomly set 1000 values to the missing value of 8
for(i in 1:1000) {
  mat8[randomRow[i], randomCol[i]] <- 8
}
stuItems8 <- mat8[,1:13]
stuItems8$oppID <- factor(rownames(mat8), levels=rownames(mat8))
stuItems8 <- reshape(data=stuItems8, varying=c(paramTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=paramTab$ItemID, timevar="key")
# it seems AM uses a missing value of c, so use that, note: weighted
mml1A <- mml(stuItems=split(stuItems8, stuItems8$oppID), stuDat=stuDat, paramTab=paramTab, Q=34, idVar="oppID", missingValue="c", weightvar="origwt")
mml1As <- summary(mml1A, varType="Taylor", stratavar="repgrp1", psuvar="jkunit", gradientHessian=TRUE)
mml1AsT_coef_REF <- structure(c(276.603996367644, 33.27212524293007,
                                0.9749103126644552, 0.9559105270104056),
                              .Dim = c(2L,2L), .Dimnames = list(c("(Intercept)", "Population SD"),
                                                                c("Estimate", "StdErr")))

expect_equal(mml1As$coefficients[,1], mml1AsT_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*20)
expect_equal(mml1As$coefficients[,2], mml1AsT_coef_REF[,2], tolerance=10*(.Machine$double.eps)^0.25)

context("missing data rows")
mat2 <- mat[1:1000,]
mat2[,1:13] <- NA
mat2 <- rbind(mat, mat2[1:1000,])
stuItems3 <- mat2[,1:13]
stuItems3$oppID <- factor(rownames(mat2), levels=rownames(mat2))
stuItems3 <- reshape(data=stuItems3, varying=c(paramTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=paramTab$ItemID, timevar="key")

rownames(mat) <- paste0("pseudo-student",1:nrow(mat))
stuDat3 <- rbind(stuDat, stuDat[1:1000,])
stuDat3$oppID <- rownames(mat2)
mml1C <- mml(stuItems=split(stuItems3, stuItems3$oppID), stuDat=stuDat3, paramTab=paramTab, Q=34, idVar="oppID"); summary(mml1C)

# test fit
expect_equal(coef(mml1), coef(mml1C), tolerance=sqrt(.Machine$double.eps)*20)
# test standard errors
mml1Cs <- summary(mml1C, gradientHessian=TRUE)
expect_equal(mml1s$coefficients, mml1Cs$coefficients, tolerance=(.Machine$double.eps)^0.25)
mml1CRobust <- summary(mml1C, varType="robust")
expect_equal(mml1CRobust$coefficients, mml1Robust$coefficients, tolerance=sqrt(.Machine$double.eps)*200)
mml1CTaylor <- summary(mml1C, varType="Taylor", stratavar="repgrp1", psuvar="jkunit")
expect_equal(mml1CTaylor$coefficients, mml1Taylor$coefficients, tolerance=sqrt(.Machine$double.eps)*20)
mml1CCluster <- summary(mml1C, varType="cluster", clustervar="repgrp1")
expect_equal(mml1CCluster$coefficients, mml1Cluster$coefficients, tolerance=sqrt(.Machine$double.eps)*20)

context("missing data rows, with regressor")

mml2C <- mml(~x1,stuItems=split(stuItems3, stuItems3$oppID), stuDat=stuDat3, paramTab=paramTab, Q=34, idVar="oppID"); summary(mml1C)

# test fit
expect_equal(coef(mml2), coef(mml2C), tolerance=sqrt(.Machine$double.eps)*20)
# test standard errors
mml2Cs <- summary(mml2C, gradientHessian=TRUE)
expect_equal(mml2s$coefficients, mml2Cs$coefficients, tolerance=(.Machine$double.eps)^0.25)
mml2CRobust <- summary(mml2C, varType="robust")
expect_equal(mml2CRobust$coefficients, mml2Robust$coefficients, tolerance=sqrt(.Machine$double.eps)*200)
mml2CTaylor <- summary(mml2C, varType="Taylor", stratavar="repgrp1", psuvar="jkunit", gradientHessian=TRUE)
expect_equal(mml2CTaylor$coefficients, mml2Taylor$coefficients, tolerance=(.Machine$double.eps)^0.25)
mml2CCluster <- summary(mml2C, varType="cluster", clustervar="repgrp1")
expect_equal(mml2CCluster$coefficients, mml2Cluster$coefficients, tolerance=sqrt(.Machine$double.eps)*200)

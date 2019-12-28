\dontrun{
# get NAEP Primer data
require(EdSurvey)

# data
sdf <- readNAEP(system.file("extdata/data", "M36NT2PM.dat", package = "NAEPprimer"))
cols <- c("m066401", "m093701", "m086001", "m051901", "m067801", "m046501",
          "origwt", "repgrp1", "jkunit", "dsex")
data <- getData(sdf, varnames=cols, addAttributes=TRUE,
                omittedLevels=FALSE, defaultConditions=FALSE,
                returnJKreplicates=FALSE)

# 3PL items only:
# P0 is the discrimination parameter (a),
# P1 is the item difficulty (d),
# P2 is the guessing parameter (g) 
# polytomous responses could use P3-P10 for more difficulties
paramTab <- structure(list(ItemID = c("m066401", "m093701", "m086001",
                                      "m051901", "m067801", "m046501"),
                           P0 = c(0.68, 1.22, 1.05, 1.6, 0.86, 1.03),
                           P1 = c(-0.33, 1.81, 1, 0.61, -1.61, -0.14),
                           P2 = c(0.15, 0.17, 0.22, 0.08, 0.06, 0.37),
                           P3 = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
                           P4 = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
                           P5 = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
                           P6 = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
                           P7 = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
                           P8 = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
                           P9 = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
                           P10 = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
                           ScorePoints = c(1, 1, 1, 1, 1, 1),
                           MODEL = c("3pl", "3pl", "3pl", "3pl", "3pl", "3pl")),
                      row.names = c(1L, 3L, 4L, 5L, 9L, 13L),
                      class = "data.frame", location = 277.1563, scale = 37.7297)
# scores an item as correct if it contains an asterisk and as skipped if it
# is "Omitted", "Not Reached", or "Multiple". The value NA is left as NA.
# this score function is intended to be simple not reflect typical NAEP scoring.
simpleScore <- function(col) {
  score0 <- 0+grepl("*", col, fixed=TRUE)
  score1 <- ifelse(col %in% c("Omitted", "Not Reached", "Multiple"), 8, score0)
  score2 <- ifelse(col %in% NA, NA, score1)
  return(score2)
}

# score each item in paramTab
for(name in paramTab$ItemID){
  # show score output vs input data
  print(table(sdf[,name], simpleScore(sdf[,name]), useNA="ifany"))
  # score item
  data[,name] <- simpleScore(data[,name])  
}

# make stuItems 
data$id <- 1:nrow(data)
# first make a long data.frame of the item score data
stuItems <- reshape(data=data, varying=c(paramTab$ItemID), idvar=c("id"),
                    direction="long", v.names="score", times=paramTab$ItemID,
                    timevar="key")[,c("id", "key", "score")]
# then break it up into a single data.frame per student
stuItems <- split(stuItems, "id")

# Studat is the student covariates, weights, and sampling information
# used for variance estimation
stuDat <- data[, c('origwt', 'repgrp1', 'jkunit', 'dsex', 'id')]

### MML call 
mml1 <- mml(~dsex, stuItems=stuItems, 
            stuDat=stuDat, paramTab=paramTab, 
            regType = 'regression', Q=34, idVar="id", weightvar = "origwt")

# summary, assumes the sample was drawn IID
summary(mml1)
# summary, accounts for correlation between students in the same schools
summary(mml1, varType="Taylor", stratavar="repgrp1", psuvar="jkunit")
}
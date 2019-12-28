#' Format AM dct File for Use with DE
#' @description
#' Takes an \code{AM dct} file and formats it for use with the \code{mml} method
#' as \code{paramTab}.
#' 
#' @param file a file location from which to read the \code{dct} file
#' 
#' @return a \code{data.frame} in a format suitable for use with \code{mml} as
#' a \code{paramTab}.
#' 
#' @author Paul Bailey
#' @export
parseAMdct <- function(file) {
  if(length(file) == 1 && file.exists(file)) {
    file <- readLines(file)
  } 
  z <- rep(NA, length(file))
  res <- data.frame(ItemID=z, P0=z, P1=z, P2=z, P3=z, P4=z,
                    P5=z, P6=z, P7=z, P8=z, P9=z, P10=z,
                    ScorePoints=z, MODEL=z, scale=z, test=z, subtest=z,
                    omitted=z,
                    stringsAsFactor=FALSE)
  parsing <- FALSE # wait for variables to parse
  startedParsing <- FALSE
  ei <- 1
  for(li in file) {
    li <- trimws(li)
    if(!parsing) {
      if(li == "Variables") {
        parsing <- TRUE
      }
    } else {
      # parsing
      if(nchar(li) > 0) {
        startedParsing <- TRUE
        print(li)
        elements <- strsplit(li, " ")[[1]]
        if(elements[1] %in% c("ORIGWT", "REPGRP1", "JKUNIt")) {
          return(res[1:(ei-1),])
        }
        res$ItemID[ei] <- elements[1]
        elements <- elements[-1]
        elementParts <- strsplit(elements, "=")
        res$MODEL[ei] <- getPart(elementParts, "IRM")
        if(tolower(res$MODEL[ei]) == "pcl") {
          res$MODEL[ei] <- "GPCM"
          # example AM data:
          #M073601 ontest=1 onsubtest=1 IRM=PCL ipa=0.58 ipb=1.13 scale=1.7 omitted=8 IPD1=-0.13 IPD2=-0.24 IPD3=0.38 IPD4=-0.01 
          # example paramTab entry:
          #    ItemID   P0   P1   P2   P3   P4 P5 P6 P7 P8 P9 P10 ScorePoints MODEL
          #   m073601 0.58 1.26 1.37 0.75 1.14 NA NA NA NA NA  NA           5  GPCM
          res$P0[ei] <- getPart(elementParts, "ipa")
          ipb <- as.numeric(getPart(elementParts, "ipb"))
          res$scale[ei] <- getPart(elementParts, "scale")
          res$omitted[ei] <- getPart(elementParts, "omitted")
          nxti <- 1
          nxt <- getPart(elementParts, paste0("IPD",nxti))
          while(!is.na(nxt)) {
            res[ei,paste0("P",nxti)] <- ipb - as.numeric(nxt)
            nxti <- nxti + 1
            nxt <- getPart(elementParts, paste0("IPD",nxti))
          }
          res$ScorePoints[ei] <- nxti
        }
        if(tolower(res$MODEL[ei]) == "3pl") {
          res$MODEL[ei] <- "3pl"
          # example AM data:
          #N202831 IRM=3pl ipa=0.75 ipb=-2.58 ipc=0.25 scale=1.7 ontest=1 onsubtest=1 omitted=8
          # example paramTab entry:
          #    ItemID   P0    P1   P2 P3 P4 P5 P6 P7 P8 P9 P10 ScorePoints MODEL
          #   n202831 0.75 -2.58 0.25 NA NA NA NA NA NA NA  NA           1   3pl
          res$P0[ei] <- getPart(elementParts, "ipa")
          res$P1[ei] <- getPart(elementParts, "ipb")
          res$P2[ei] <- getPart(elementParts, "ipc")
          res$scale[ei] <- getPart(elementParts, "scale")
          res$omitted[ei] <- getPart(elementParts, "omitted")
          res$ScorePoints[ei] <- 1
        }
        res$test[ei] <- getPart(elementParts, "ontest")
        res$subtest[ei] <- getPart(elementParts, "onsubtest")
        ei <- ei + 1
      } else {
        if(startedParsing) {
          parsing <- FALSE
        }
      }
    }
  }
  res <- res[1:(ei-1),]
  res$ItemID <- tolower(res$ItemID)
  return(res)
}

getPart <- function(le, element) {
  res <- lapply(le, function(e) {
    if(tolower(e[1]) == tolower(element)) {
      return(e[2])
    }
    return(NULL)
  })
  res <- unlist(res)
  if(is.null(res)) {
    return(NA)
  }
  suppressWarnings(rr <- as.numeric(res))
  if(is.na(rr)) {
    return(res)
  }
  return(rr)
}

csv2paramList <- function(pModel = c('GPCM', 'GRM', 'PCM'), dat, subsetItems = NULL){
  # If converting modelParameters tab from config, only grab
  # relevant variables and rename ItemID to ITEM_ID
  if (all(c('MODEL', 'ItemID', 'ScorePoints', paste0('P', 0:10)) %in% names(dat))) {
    dat <- dat[, c('ItemID', 'MODEL', 'ScorePoints', paste0('P', 0:10))]
    names(dat)[1] <- 'ITEM_ID'
  }
  
  params2use <- with(dat, dat[order(dat$ITEM_ID), ])

  if(!is.null(subsetItems)){
    dat <- dat[which(dat$ITEM_ID %in% subsetItems), ]
    params2use <- dat[match(subsetItems, dat$ITEM_ID), ] #sort according to item order input
  }
  
  pModel <- toupper(pModel)
  
  names(params2use) <- tolower(names(params2use))
  finalItems <- params2use$item_id
  
  # any polytymous items?-- this will take care of IRTGPC sp = 2, IRTGPC sp = 3, IRTPCL sp = 2, IRTPCL sp = 3
  any.poly <- any(params2use$scorepoints > 1)
  # identifly all other items-- IRTGPC sp = 1, IRT3pln sp = 1, IRTPCL sp = 1, IRT3pl sp = 1
  # in config file, IRTPCL scorepoints = 1 items have b parameter stored in p0, unlike other 3pl items
  # for example, parameters for IRT3pl sp = 1 are stored as p0 = a(1), p1 = b, p2 = c(0)
  pos <- which(params2use$scorepoints == 1 & grepl('3pl', tolower(params2use$model))) # 3pl
  pos1 <- which(params2use$scorepoints == 1 & tolower(params2use$model) %in% 'irtgpc') # gpc
  pos <- c(pos, pos1)
  pos2 <- which(params2use$scorepoints == 1 & tolower(params2use$model) %in% 'irtpcl') # pcl
  
  mcParams <- params2use[pos, paste0('p', 0:2)]
  mcParams$p2 <- with(mcParams, ifelse(is.na(p2), 0, p2)) # make c-parameter if it is missing, this is needed for IRTGPC sp = 1
  
  mcParams2 <- params2use[pos2, paste0('p', 0:2)]
  if(nrow(mcParams2) > 0){
    mcParams2$p1 <- mcParams2$p0 # b is stored in p0 for IRTPCL sp = 1 items
    mcParams2$p0 <- 1 # set a to 1
    mcParams2$p2 <- 0 # set c to 0
  }
  
  mcParams <- rbind(mcParams, mcParams2)
  binary <- c(pos, pos2)
  
  ## inserted 4 lines for modifying
  binaryID <- params2use$item_id[binary]  
  polyID <- params2use$item_id[-binary]
  
  if(any.poly){
    ### This is the code to get the CR items organized for the params list
    numSteps <- length( grep("[p][0-9]", colnames(params2use)) )
    if(pModel == 'PCM') steps <- paste('p', 0:(numSteps-1), sep='')
    if(pModel %in% c('GRM', 'GPCM')) steps <- paste('p', 1:(numSteps-1), sep='')
    if(length(binary) > 0){
      crParams <- params2use[-binary, steps]
    } else {
      crParams <- params2use[, steps]
    } 
    if(pModel == 'PCM') aVar <- rep(1, nrow(crParams))
    if(length(binary) > 0){
      if(pModel %in% c('GRM', 'GPCM')) aVar <- params2use[-binary, 'p0']
    } else{
      if(pModel %in% c('GRM', 'GPCM')) aVar <- params2use[, 'p0']
    }
    crParams$base <- 0
    if(pModel == 'PCM') crParams <- crParams[, c('base', steps)]
    if(pModel == 'GPCM') crParams <- crParams[, c('base', steps)]
    if(pModel == 'GRM') crParams <- crParams[, steps]
    crParams <- as.data.frame(t(crParams))
    crList <- as.list(crParams)
    crList <- lapply(crList, function(x) x[!is.na(x)])
    if(pModel == 'PCM') pModel <- 'GPCM' # so that model can be automatically called by irt.ability
    if(length(binary) > 0){ # test has both CM and polytomous items
      params <- list('3pl' = list(a = mcParams$p0, b = mcParams$p1, c = mcParams$p2), 
                     gpcm = list(a = aVar , d = crList, polyID=polyID), model = pModel, items = finalItems, mcID=binaryID, pyID=polyID, mcPos = binary, 
                     score.pts = params2use$scorepoints)
    } else { # test has only polytomous items
      params <- list('3pl' = NULL, 
                     gpcm = list(a = aVar , d = crList), model = pModel, items = finalItems, mcID=NULL, pyID=polyID, mcPos = NULL, score.pts = params2use$scorepoints)
    }
  } else { # test has only MC items
    if(pModel == 'PCM') pModel <- 'GPCM' # so that model can be automatically called by irt.ability
    params <- list('3pl' = list(a = mcParams$p0, b = mcParams$p1, c = mcParams$p2), 
                   gpcm = NULL, items = finalItems, mcID=binaryID, pyID=NULL, mcPos = binary, score.pts = params2use$scorepoints, model = pModel)
  }
  class(params) <- 'param.list'
  params     
}

getPos <- function(...) UseMethod("getPos") 
getPos.param.list <- function(object, ...){
  object$mcPos
}

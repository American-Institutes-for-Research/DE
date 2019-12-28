#' @importFrom foreach %dopar% foreach
#' @importFrom iterators iter
calcRR1_dopar <- function(stu, Q, model, paramTab, missingCode, missingValue, D, nodes) {
  itx1 <- iter(stu)
  rr1 <- foreach(dopari = itx1, .export = c('csv2paramList', 'getPos', 'getPos.param.list', 'ldbinom2', "gpcm", "grm")) %dopar% {
    x <- dopari$score
    if(all(is.na(x))) {
      return(rep(0,Q))
    }
    reqVar <- c('MODEL', 'ItemID', 'ScorePoints', paste0('P', 0:10))
    missIndex <- which(!reqVar %in% names(paramTab))
    # check for missingness
    if (length(missIndex) > 0) {
      stop(paste('Missing variable(s):', reqVar[missIndex]))
    }
    params <- csv2paramList(toupper(model), paramTab, dopari$key)
    # get parameter values from items
    a <- params$'3pl'$a
    b <- params$'3pl'$b
    c <- params$'3pl'$c
    aa <- params$gpcm$a
    d <- params$gpcm$d
    ind.dichot <- getPos(params)
    x1 <- x[ind.dichot]
    if(missingValue %in% "c") {
      x1[x1 %in% missingCode] <- c[x1 %in% missingCode] # values of missingCode get recoded to the guessing parameter
    } else {
      x1[x1 %in% missingCode] <- missingValue # values of missingCode get recoded to "missingValue" correct
    }
    # fix multinomial items
    if (model == 'gpcm') {
      x2 <- x[!1:length(x) %in% ind.dichot] +1 
    } else {
      x2 <- x[!1:length(x) %in% ind.dichot]
    }
    # get multiple choice part
    log.mc.part <- if (length(x1) > 0) {
      sapply(1:Q, function(j) sum(ldbinom2(x1, c + (1 - c) / (1 + exp(-D * a * (nodes[j] - b)))), na.rm=TRUE))
    } else {
      0
    }
    # get polytimous part
    log.poly.part <- if (length(x2) > 0) {
      sapply(1:Q, function(j) sum(log(mapply(get(model), d = d, a = aa, theta = nodes[j], score = x2, D=D)), na.rm=TRUE))
    } else {
      0
    }
    log.parts <- rbind(log.mc.part, log.poly.part)
    return(exp(apply(log.parts, 2, sum, na.rm=TRUE)))
  }
  return(do.call(cbind, rr1))
}

calcRR1 <- function(stu, Q, model, paramTab, missingCode, missingValue, D, nodes) {
  rr1 <- lapply(stu, function(dopari) {
    x <- dopari$score
    if(all(is.na(x))) {
      return(rep(0,Q))
    }
    reqVar <- c('MODEL', 'ItemID', 'ScorePoints', paste0('P', 0:10))
    missIndex <- which(!reqVar %in% names(paramTab))
    # check for missingness
    if (length(missIndex) > 0) {
      stop(paste('Missing variable(s):', reqVar[missIndex]))
    }
    params <- csv2paramList(toupper(model), paramTab, dopari$key)
    # get parameter values from items
    a <- params$'3pl'$a
    b <- params$'3pl'$b
    c <- params$'3pl'$c
    aa <- params$gpcm$a
    d <- params$gpcm$d
    ind.dichot <- getPos(params)
    x1 <- x[ind.dichot]
    if(missingValue %in% "c") {
      x1[x1 %in% missingCode] <- c[x1 %in% missingCode] # values of missingCode get recoded to the guessing parameter
    } else {
      x1[x1 %in% missingCode] <- missingValue # values of missingCode get recoded to "missingValue" correct
    }
    # fix multinomial items
    if (model == 'gpcm') {
      x2 <- x[!1:length(x) %in% ind.dichot] +1 
    } else {
      x2 <- x[!1:length(x) %in% ind.dichot]
    }
    # get multiple choice part
    log.mc.part <- if (length(x1) > 0) {
      sapply(1:Q, function(j) sum(ldbinom2(x1, c + (1 - c) / (1 + exp(-D * a * (nodes[j] - b)))), na.rm=TRUE))
    } else {
      0
    }
    # get polytimous part
    log.poly.part <- if (length(x2) > 0) {
      sapply(1:Q, function(j) sum(log(mapply(get(model), d = d, a = aa, theta = nodes[j], score = x2, D=D)), na.rm=TRUE))
    } else {
      0
    }
    log.parts <- rbind(log.mc.part, log.poly.part)
    return(exp(apply(log.parts, 2, sum, na.rm=TRUE)))
  })
  return(unname(do.call(cbind, rr1)))
}

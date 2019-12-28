#' Marginal Maximum Likelihood Estimation of Linear Models
#' @description
#' Implements a survey-weighted marginal maximum estimation, a type of
#' regression where the outcome is a latent trait (such as student ability.
#' Instead of using an estimate, the likelihood function marginalizes student
#' ability. Includes a variety of variance estimation strategies.
#' 
#' @param formula  a formula object in the style of \code{lm}
#' @param stuItems a list where each element is named a student ID and contains
#'                 a \code{data.frame}; see Details for the format
#' @param stuDat a \code{data.frame} with a single row per student. Predictors in
#'               the formula must be in \code{stuDat}.
#' @param paramTab a \code{data.frame} with columns shown in Details
#' @param Q the number of integration points
#' @param polyModel polytomous response model;
#'                  one of \code{GPCM} for the Graded Partial Credit Model
#'                  or \code{GRM} for the Graded Response Model
#' @param regType one of \code{regression} or \code{popMean} where the latter
#'                estimates a population level mean
#' @param weightvar a variable name on \code{stuDat} that is the full sample weight
#' @param control a list with four elements that control the fitting process. See Details.
#' @param idVar a variable name on \code{stuDat} that is the identifier. Every 
#'              ID from \code{stuDat} must appear on \code{stuItems} and vice versa.
#' @param missingCode the value a score is set to that indicates the item is missing.
#'                    An item scored as \code{NA} will be ignored. The \code{missingCode}
#'                    argument allows the user to recode scores to \code{missingValue}.
#'                    This argument applies exclusively to binomial items.
#' @param missingValue the value to set items scored as \code{missingCode}.
#'                     When set to a number, that value is used for all items.
#'                     When set to \dQuote{\code{C}}, then the guessing parameter
#'                     is used.
#' @param multiCore allows the \code{foreach} package to be used. You should
#'                  have already called 
#' \ifelse{latex}{the \code{registerDoParallel} function in the \code{doParallel} package}{\code{\link[doParallel]{registerDoParallel}}}.
#' @param bobyqaControl a list that gets passed to \ifelse{latex}{the \code{bobyqa} optimizer in \code{minqa}}{\code{\link[minqa]{bobyqa}}}
#' @details
#' 
#' The \code{mml} function models a latent outcome conditioning on student
#' item response data, student covariate data, and item parameter information;
#' these three parts are broken up into three arguments.
#' Student item response data go into \code{stuItems}, whereas student
#' covariates, weights, and sampling information go into \code{stuDat}.
#' The \code{paramTab}
#' contains item parameter information for each item---the result of a
#' separate item parameter scaling. In the case of 
#' the National Assessment of Educational Progress (NAEP),
#' they can be found online, for example, at
#' \href{https://nces.ed.gov/nationsreportcard/tdw/analysis/scaling_irt.aspx}{https://nces.ed.gov/nationsreportcard/tdw/analysis/scaling_irt.aspx}.
#' The model for dichotomous responses data is by default three Parameter Logit
#' (3PL), unless the item parameter information provided by users suggests
#' otherwise. For example, if the scaling used a two Parameter Logit (2PL) model,
#' then the guessing parameter can simply be set to zero. For polytomous
#' responses data, the model is dictated by the \code{polyModel} argument.
#' 
#' Student data are broken up into two parts. The item response data goes
#' into \code{stuItems} ,and the student covariates for the formula go into
#' \code{stuDat}. Information about items, such as item difficulties, is in 
#' \code{paramTab}. All dichotomous items are assumed to be 
#' 3PL, though by setting the guessing parameter to zero, the user
#' can use a 2PL or the one Parameter Logit (1PL) or Rasch models.
#' The model for polytomous responses data is dictated by the \code{polyModel}
#' argument.
#' 
#' The marginal maximum likelihood then integrates the product of the student
#' ability from the assessment data, and the estimate from the linear model
#' estimates each student's ability based on the \code{formula} provided
#' and a residual standard error term. This integration happens from the
#' minimum node to the maximum node in the \code{control} argument (described
#' later in this section) with \code{Q} quadrature points. 
#' 
#' The \code{stuItems} argument has the scored student data. It is a list where
#' each element is named with student ID and contains
#' a \code{data.frame} with at least two columns.
#' The first required column is named
#' \code{key} and shows the item name as it appears in \code{paramTab};
#' the second column in named
#' \code{score} and shows the score for that item. For binomial
#' items, the \code{score} is 0 or 1. For \code{GPCM} items, the scores
#' start at zero as well. For \code{GRM}, the scores start at 1.
#' 
#' The \code{paramTab} argument is a \code{data.frame} with a column named
#' \code{ItemID} that agrees with
#' the \code{key} column in the \code{stuItems} argument,
#' and, for  a 3PL item, columns \code{P0},
#' \code{P1}, and \code{P2} for the \dQuote{a}, \dQuote{d}, and
#' \dQuote{g} parameters, respectively; see the vignette for details of
#' the 3PL model.
#' For a \code{GPCM} model, \code{P0} is the \dQuote{a} parameter, and the other 
#' columns are the \dQuote{d} parameters; see the vignette for details
#' of the GPCM model.
#' 
#' The \code{control} argument is a list with, optional, items \code{D}, the
#' scale parameter, that defaults to 1.7; \code{startVal}, which is the starting
#' value for the coefficients; and \code{min.node} and \code{max.node}, which
#' sets the range of nodes for all students; these default to
#' -4 and 4, respectively. The quadrature points then are a range
#' from \code{min.node} to \code{max.node} with a total of \code{Q} nodes.
#' 
#' @return object of class \code{mml.means}. 
#' This is a list with elements: 
#' \item{call}{the call used to generate this \code{mml.means} object}
#' \item{coefficients}{the marginal maximum likelihood regression coefficients, including the estimated residual standard error}
#' \item{LogLik}{the log-likelihood of the fit model}
#' \item{X}{the design matrix of the marginal maximum likelihood regression}
#' \item{Convergence}{a convergence note from the \code{bobyqa} optimizer}
#' \item{location}{used for scaling the estimates}
#' \item{scale}{used for scaling the estimates}
#' \item{lnlf}{the likelihood function} 
#' \item{rr1}{the density function of each individual, conditional only on item responses in \code{stuItems}}
#' \item{stuDat}{the \code{stuDat} argument}
#' \item{weightvar}{the weight variable}
#' \item{nodes}{the nodes the likelihood was evaluated on}
#' \item{iterations}{the number of iterations required to reach convergence}
#' \item{obs}{the number of observations used}
#' 
#' @example \man\examples\de.R
#' @author Harold Doran, Paul Bailey, Claire Kelley, and Sun-joo Lee 
#' @export
#' @importFrom minqa bobyqa 
#' @importFrom stats as.formula model.matrix coef dbinom sd
#' @importFrom utils head
mml <- function(formula, stuItems, stuDat, paramTab, Q=30, polyModel=c('GPCM', 'GRM'),
                regType=c('regression', 'popMean'), 
                weightvar=NULL,
                control=list(),
                idVar=c(),
                missingCode=8,
                missingValue="c",
                multiCore=FALSE,
                bobyqaControl=list) { # missing code from NAEP
  call <- match.call()
  polyModel <- match.arg(polyModel)
  regType <- match.arg(regType)
  polyModel <- tolower(polyModel)
  if(!idVar %in% names(stuDat)) {
    stop(paste0(dQuote("idVar"), " must be a variable on stuDat to confirm the ", dQuote("stuDat"), " rows agree with the ", dQuote("stuItems"), " list."))
  }
  if(!inherits(stuDat, "data.frame")) {
    stop(paste0("Argument ", dQuote("stuDat"), " must be a data frame."))
  } else {
    if(length(class(stuDat))>1) {
      # not all data.frames are the same, recast
      stuDat <- as.data.frame(stuDat)
    }
  }
  if(!inherits(paramTab, "data.frame")) {
    stop(paste0("Argument ", dQuote("paramTab"), " must be a data frame."))
  } else {
    if(length(class(paramTab))>1) {
      # not all data.frames are the same, recast
      paramTab <- as.data.frame(paramTab)
    }
  }
  regType <- match.arg(regType)
  if(any(!names(stuItems) %in% stuDat[[idVar]])) {
    missing <- names(stuItems)[!names(stuItems) %in% stuDat[[idVar]]]
    stop(paste0("The ", dQuote("stuItems"), " argument must be a list with names that correspond to every " , dQuote("idVar"), " in ", dQuote("stuDat"), ". some missing IDs ", pasteItems(dQuote(head(missing,5))), "."))
  }
  if(any(!stuDat[[idVar]] %in% names(stuItems))) {
    missing <- stuDat[[idVar]][!stuDat[[idVar]] %in% names(stuItems)]
    stop(paste0("The ", dQuote("stuDat"), " argument must be a data.frame with a column ", dQuote("idVar"), " that correspond to every name of ", dQuote("stuItems"), ". some missing IDs ", pasteItems(dQuote(head(missing,5))), "." ))
  }
  # make sure stuItems and stuDat are in the same order
  stuItems <- stuItems[stuDat[[idVar]]]
  
  # make sure stuItems (now stu) has only data.frames in each element
  stu <- lapply(stuItems, function(x) {
    return(as.data.frame(x)[c(idVar, 'key', 'score')])
  })

  # subset to students who have at least one valid score
  stuValid <- unlist(lapply(stu, function(x) {
    return(!all(is.na(x$score)))
  }))
  stu <- stu[stuValid]
  stuDat <- stuDat[stuValid,]

  # now form X
  if (regType == 'regression') {
    if(missing(formula)) {
      formula <- ~ 1
    } # end if(missing(formula))
    X <- model.matrix(formula, stuDat)
    K <- ncol(X) # number of fixed parameters to estimate   
    nms <- c(colnames(X), 'Population SD') 
    startVal <- c(rep(0, K), 1)
    if(nrow(X) != nrow(stuDat)) {
      stop("Missing not allowed in independent variables.")
    }
  }
  if (regType == 'popMean') {
    startVal <- c(0, 1)
    nms <- c('Population Mean', 'Population SD')
  }   

  # setup controls, including startVal
  con <- list(D = 1.7, startVal = startVal, min.node = -4, max.node = 4)
  con[names(control)] <- control
  ### These are the exported functions passed to the dopar function
  nodes <- seq(from = con$min.node, to = con$max.node, length.out = Q)

  ### This portion of the code computes all the likelihood evaluations
  ### and does so outside of the function that is maximized
  ### This saves a lot of overhead by using fixed quadrature points
  if(multiCore) {
    rr1 <- calcRR1_dopar(stu, Q, polyModel, paramTab, missingCode, missingValue, con$D, nodes)
  } else {
    rr1 <- calcRR1(stu, Q, polyModel, paramTab, missingCode, missingValue, con$D, nodes)
  }

  if (regType == 'popMean') {
    regType <- "regression"
    X <- matrix(1, nrow=ncol(rr1), ncol=1)
  }
  if (regType == 'regression') {
    fn2 <- fn.regression(X_=X, i=NULL, wv=weightvar, rr1=rr1, nodes=nodes, stuDat=stuDat)
    opt <- bobyqa(con$startVal, fn2)
    opt$par[length(opt$par)] <- abs(opt$par[length(opt$par)])
    iter <- opt$feval
    convergence <- opt$msg
    names(opt$par) <- c(colnames(X), "s")
  }
  # get location and scale from paramTab
  location <- 0
  scale <- 1
  if(!is.null(attr(paramTab, "location"))) {
    location <- attr(paramTab, "location")
    if(!is.null(attr(paramTab, "scale"))) {
      scale <- attr(paramTab, "scale")
    } else {
      scale <- 1
    }
  } else {
    location <- 0
    scale <- 1
  }
  if(is.null(weightvar)) {
    obs <- sum(apply(rr1,2,sum)>0)
  } else {
    obs <- sum(stuDat[,weightvar]>0 & apply(rr1,2,sum)>0) # number with positive weight
  }
  coefficients <- opt$par
  names(coefficients) <-nms
  structure(list("call" = call,
                 "coefficients" = coefficients,
                 "LogLik" = -1/2*opt$fval,
                 "X" = X,
                 'Convergence' = convergence,
                 "location" = location,
                 "scale" = scale,
                 "lnlf"= fn2,
                 "rr1"= rr1,
                 "stuDat" = stuDat,
                 "weightvar" = weightvar,
                 "nodes" = nodes,
                 iterations = iter,
                 "obs" = obs),
            class = "mml.means")
}

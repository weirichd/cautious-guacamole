#' fits a 2-dimensional thin-plate spline smoothing spline anova model to the elements of the
#' cholesky factor when the innovation variances are known.
#'
#' @param formula
#' @param type List specifying the type of spline for each variable.
#' @param data list containing the data where the matrix of 2-dimensional regressors are a single element in the list
#' @param weights vector of the N x M-1 innovation variances
#' @param subset Optional vector specifying a subset of observations to be used in the fitting process.
#' @param offset Optional offset term with known parameter 1.
#' @param na.action Function which indicates what should happen when the data contain NAs.
#' @param partial Optional symbolic description of parametric terms in partial spline models.
#' @param method Method for smoothing parameter selection. Supported are method="v" for GCV, method="m" for GML (REML), and method="u" for Mallows' CL.
#' @param alpha	Parameter modifying GCV or Mallows' CL; larger absolute values yield smoother fits; negative value invokes a stable and more accurate GCV/CL evaluation algorithm but may take two to five times as long. Ignored when method="m" are specified.
#' @param varht external variance estimate needed for method="u". Ignored when method="v" or method="m" are specified.
#' @param id.basis Index designating selected "knots".
#' @param nbasis number of "knots" to be selected. Ignored when id.basis is supplied.
#' @param seed Seed to be used for the random generation of "knots". Ignored when id.basis is supplied.
#' @param random	Input for parametric random effects in nonparametric mixed-effect models. See mkran for details.
#' @param skip.iter Flag indicating whether to use initial values of theta and skip theta iteration. See notes on skipping theta iteration.
#' @return a list of class vcssanova
#' @seealso
#' @export
#' @examples

vcssanova <- function (formula=as.formula("y~x"),
                       type = NULL, wt, subset=NULL,
                       data, offset=NULL, na.action = na.omit,
                       partial = NULL, method = "v",
                    alpha = 1.4, varht = 1,
                    nbasis = NULL, seed = NULL,
                    random = NULL, skip.iter = FALSE)
{
  mf <- match.call()
  mf$type <- mf$method <- mf$varht <- mf$partial <- NULL
  mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
  mf$random <- mf$skip.iter <- NULL

  data$x <- as.matrix(data$x)
  dimnames(data$x)[1:2] <- NULL
  mf$formula <- as.formula("y~x")
  mfr <- model.frame(formula=mf$formula,
                     data=list(x=data$x,
                               y=rep(1,dim(data$x)[1])))
  mf <- mfr
  rm(mfr)
  nobs <- dim(mf)[1]
  id.basis <- 1:nobs
  if (is.null(id.basis)) {
    if (is.null(nbasis))
      nbasis <- max(30, ceiling(10 * nobs^(2/9)))
    if (nbasis >= nobs)
      nbasis <- nobs
    if (!is.null(seed))
      set.seed(seed)
    id.basis <- sample(nobs, nbasis, prob = wt)
  }
  else {
    if (max(id.basis) > nobs | min(id.basis) < 1)
      stop("gss error in ssanova: id.basis out of range")
    nbasis <- length(id.basis)
  }
  term <- mkterm(mf, type)
    if (!is.null(random)) {
    if (class(random) == "formula")
      random <- mkran(random, data)
  }
  s <- q <- NULL
  nq <- 0
  for (label in term$labels) {
    if (label == "1") {
      s <- cbind(s, rep(1, len = nobs))
      next
    }
    x <- mf[, term[[label]]$vlist]
    x.basis <- mf[id.basis, term[[label]]$vlist]
    nphi <- term[[label]]$nphi
    nrk <- term[[label]]$nrk
    if (nphi) {
      phi <- term[[label]]$phi
      for (i in 1:nphi) s <- cbind(s, phi$fun(x, nu = i,
                                              env = phi$env))
    }
    if (nrk) {
      rk <- term[[label]]$rk
      for (i in 1:nrk) {
        nq <- nq + 1
        q <- array(c(q, rk$fun(x, x.basis, nu = i, env = rk$env,
                               out = TRUE)), c(nobs, nbasis, nq))
      }
    }
  }
  if (is.null(q)) {
        stop("gss error in ssanova: use lm for models with only unpenalized terms")        
  }
  if (!is.null(partial)) {
    mf.p <- model.frame(partial, data)
    for (lab in colnames(mf.p)) mf[, lab] <- mf.p[, lab]
    mt.p <- attr(mf.p, "terms")
    lab.p <- labels(mt.p)
    matx.p <- model.matrix(mt.p, data)[, -1, drop = FALSE]
    if (dim(matx.p)[1] != dim(mf)[1])
      stop("gss error in ssanova: partial data are of wrong size")
    matx.p <- scale(matx.p)
    center.p <- attr(matx.p, "scaled:center")
    scale.p <- attr(matx.p, "scaled:scale")
    s <- cbind(s, matx.p)
    part <- list(mt = mt.p, center = center.p, scale = scale.p)
  }
  else part <- lab.p <- NULL
  if (qr(s)$rank < dim(s)[2]){
        stop("gss error in ssanova: unpenalized terms are linearly dependent")        
  }

  W <- matrix(data=0,nrow=nrow(data$y)*(ncol(data$y)-1),
              ncol=choose(ncol(data$y),2))
  no.skip <- 0
  for (t in 2:ncol(data$y)){
        W[((0:(nrow(data$y)-1))*(ncol(data$y)-1)) + t-1,(no.skip+1):(no.skip+t-1)] <- data$y[,1:(t-1)]
        no.skip <- no.skip + t - 1
  }

  W <- matrix(data=0,nrow=nrow(data$y)*(ncol(data$y)-1),
              ncol=choose(ncol(data$y),2))
  no.skip <- 0
  for (t in 2:ncol(data$y)){
        W[((0:(nrow(data$y)-1))*(ncol(data$y)-1)) + t-1,
          (no.skip+1):(no.skip+t-1)] <- data$y[,1:(t-1)]
        no.skip <- no.skip + t - 1
  }
  
  y <- as.vector(t(data$y[,-1]))
  Dinv <- diag(1/wt)
  
  if (!is.null(offset)) {
        term$labels <- c(term$labels, "offset")
        term$offset <- list(nphi = 0, nrk = 0)
        y <- y - offset
  }
  ## ------------------------------------------------------
  M <- t(W) %*% Dinv %*% W
  Minv <- solve(M)
  y <- t(W) %*% Dinv %*% y
  if(nq==1){
        q[,,1] <- M %*% as.matrix(q[1:dim(q)[1],1:dim(q)[2],1]) %*% M
  }
  s <- M %*% s
  if (!is.null(offset)) {
    term$labels <- c(term$labels, "offset")
    term$offset <- list(nphi = 0, nrk = 0)
    y <- y - offset
  }
  ## Fit the model
  if (nq==1) {
        q <- q[,,1]
        z <- vcsspreg0(s,q,y,M,method,varht)
  }
  else z <- mspreg0(s,q,y,method,varht,prec,maxiter)
  ## Brief description of model terms
  desc <- NULL
  for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
  if (!is.null(partial)) {
        desc <- rbind(desc,matrix(c(1,0),length(lab.p),2,byrow=TRUE))
  }
  desc <- rbind(desc,apply(desc,2,sum))
  if (is.null(partial)) rownames(desc) <- c(term$labels,"total")
  else rownames(desc) <- c(term$labels,lab.p,"total")
  colnames(desc) <- c("Unpenalized","Penalized")
  ## Return the results
  obj <- c(list(call=match.call(),mf=mf,terms=term,partial=part,lab.p=lab.p,
                desc=desc),z)
  class(obj) <- c("vcssanova0","vcssanova")
  obj
}


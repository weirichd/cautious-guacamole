get_vcssanova_basis <- function (formula=as.formula("y~x"), data,
                                 type = NULL, wt, subset=NULL,
                                 offset=NULL, na.action = na.omit,
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
#      list(formula=mf$formula,x=data$x,y=rep(1,dim(data$x)[1]),mf=mf)
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
      list(y=y,M=M,q=q,s=s,nq=nq,q=q)
}
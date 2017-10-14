
## Fit Single Smoothing Parameter REGression
vcsspreg0 <- function(s,q,y,M, method="v",varht=1)
{
  ## Check inputs
  if (is.vector(s)) s <- as.matrix(s)
  if (!(is.matrix(s)&is.matrix(q)&is.vector(y)&is.character(method))) {
    stop("gss error in sspreg: inputs are of wrong types")
  }
  nobs <- length(y)
  nnull <- dim(s)[2]
  if (!((dim(s)[1]==nobs)&(dim(q)[1]==nobs)&(dim(q)[2]==nobs)
        &(nobs>=nnull)&(nnull>0))) {
    stop("gss error in sspreg: inputs have wrong dimensions")
  }
  ## Set method for smoothing parameter selection
  code <- (1:3)[c("v","m","u")==method]
  if (!length(code)) {
    stop("gss error: unsupported method for smoothing parameter selection")
  }
  # vmu, s, lds, nobs, nnull, y, z, M, q,
  # *ldq, tol, job, limnla, nlaht, score, varht, c, d, qraux, jpvt,
  # *wk, twk, twk2, qwk, info
  z <- .Fortran("dsidr0",
                as.integer(code),
                swk=as.double(s), as.integer(nobs),
                as.integer(nobs), as.integer(nnull),
                as.double(y), double(nobs), as.double(M),
                qwk=as.double(q), as.integer(nobs),
                as.double(0), as.integer(0), double(2),
                nlambda=double(1), score=double(1), varht=as.double(varht),
                c=double(nobs), d=double(nnull),
                qraux=double(nnull), jpvt=integer(nnull),
                double(2*nobs),double(2*nobs), double(2*nobs),
                as,double(q),
                info=integer(1),PACKAGE="gss")
  ## Check info for error
  if (info<-z$info) {               
    if (info>0)
      stop("gss error in sspreg: matrix s is rank deficient")
    if (info==-2)
      stop("gss error in sspreg: matrix q is indefinite")
    if (info==-1)
      stop("gss error in sspreg: input data have wrong dimensions")
    if (info==-3)
      stop("gss error in sspreg: unknown method for smoothing parameter selection.")
  }
  ## Return the fit
  c(list(method=method,theta=0),
    z[c("c","d","nlambda","score","varht","swk","qraux","jpvt","qwk")])
}

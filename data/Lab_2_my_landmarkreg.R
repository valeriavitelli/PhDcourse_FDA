mylandmarkreg <- function (unregfd, ximarks, x0marks, x0lim = NULL, WfdPar = NULL, 
          WfdPar0 = NULL, ylambda = 1e-10) 
{
  if (!(inherits(unregfd, "fd"))) 
    stop("Argument unregfd  not a functional data object.")
  Ybasis <- unregfd$basis
  nbasis <- Ybasis$nbasis
  rangeval <- Ybasis$rangeval
  if (is.null(x0lim)) 
    x0lim = rangeval
  if (is.numeric(ximarks)) {
    nximarks <- length(ximarks)
    if (is.vector(ximarks)) 
      ximarks <- matrix(ximarks, 1, nximarks)
    if (is.data.frame(ximarks)) 
      ximarks <- as.matrix(ximarks)
  }
  else {
    stop("Argument ximarks is not numeric.")
  }
  if (is.numeric(x0marks)) {
    nx0marks <- length(x0marks)
    if (is.vector(x0marks)) 
      x0marks <- matrix(x0marks, 1, nx0marks)
  }
  else {
    stop("Argument x0marks is not numeric.")
  }
  if (ncol(ximarks) != length(x0marks)) 
    stop("The number of columns in ximarks is not equal to length of x0marks.")
  if (any(ximarks <= rangeval[1]) || any(ximarks >= rangeval[2])) 
    stop("Argument ximarks has values outside of range of unregfd.")
  if (any(x0marks <= x0lim[1]) || any(x0marks >= x0lim[2])) 
    stop("Argument x0marks has values outside of range of target interval.")
  ncurve <- dim(ximarks)[1]
  if (is.null(WfdPar)) {
    Wnbasis <- length(x0marks) + 2
    Wbasis <- create.bspline.basis(rangeval, Wnbasis)
    Wfd <- fd(matrix(0, Wnbasis, 1), Wbasis)
    WfdPar <- fdPar(Wfd, 2, 1e-10)
  }
  else {
    WfdPar <- fdParcheck(WfdPar, 1)
    Wfd <- WfdPar$fd
    Wbasis <- Wfd$basis
    Wnbasis <- Wbasis$nbasis
  }
  if (is.null(WfdPar0)) {
    Wnbasis0 <- length(x0marks) + 3
    Wbasis0 <- create.bspline.basis(x0lim, Wnbasis0)
    Wfd0 <- fd(matrix(0, Wnbasis0, 1), Wbasis0)
    WfdPar0 <- fdPar(Wfd0, 2, 1e-10)
  }
  else {
    WfdPar0 <- fdParcheck(WfdPar0, 1)
    Wfd0 <- WfdPar0$fd
    Wbasis0 <- Wfd0$basis
    Wnbasis0 <- Wbasis0$nbasis
  }
  nfine <- min(c(101, 10 * nbasis))
  xfine <- seq(rangeval[1], rangeval[2], length = nfine)
  xfine0 <- seq(x0lim[1], x0lim[2], length = nfine)
  yfine <- eval.fd(xfine, unregfd)
  yregmat <- yfine
  hfunmat <- matrix(0, nfine, ncurve)
  hinvmat <- matrix(0, nfine, ncurve)
  xval <- matrix(c(x0lim[1], x0marks, x0lim[2]), nx0marks + 
                   2, 1)
  Wcoef <- matrix(0, Wnbasis, ncurve)
  nval <- length(xval)
  if (ncurve > 1) 
    cat("Progress:  Each dot is a curve\n")
  for (icurve in 1:ncurve) {
    if (ncurve > 1) 
      cat(".")
    yval <- matrix(c(rangeval[1], ximarks[icurve, ], rangeval[2]), 
                   nx0marks + 2, 1)
    Wfd <- smooth.morph(xval, yval, rangeval, WfdPar)$Wfdobj
    hfun <- monfn(xfine, Wfd)
    b <- (rangeval[2] - rangeval[1])/(hfun[nfine] - hfun[1])
    a <- rangeval[1] - b * hfun[1]
    hfun <- a + b * hfun
    hfun[c(1, nfine)] <- rangeval
    Wcoefi <- Wfd$coef
    Wcoef[, icurve] <- Wcoefi
    hfunmat[, icurve] <- hfun
    Wcoefi <- Wfd$coefs
    Wfdinv <- smooth.morph(hfun, xfine, x0lim, WfdPar0)$Wfdobj
    hinv <- monfn(xfine, Wfdinv)
    b <- (x0lim[2] - x0lim[1])/(hinv[nfine] - hinv[1])
    a <- x0lim[1] - b * hinv[1]
    hinv <- a + b * hinv
    hinv[c(1, nfine)] <- rangeval
    hinvmat[, icurve] <- hinv
    yregfd <- smooth.basis(hinv, yfine[, icurve], Ybasis)$fd
    yregmat[, icurve] <- eval.fd(xfine, yregfd, 0)
  }
  if (ncurve > 1) 
    cat("\n")
  regfdPar <- fdPar(Ybasis, 2, ylambda)
  regfd <- smooth.basis(xfine, yregmat, regfdPar)$fd
  regnames <- unregfd$fdnames
  names(regnames)[3] <- paste("Registered", names(regnames)[3])
  regfd$fdnames <- regnames
  warpfd <- smooth.basis(xfine, hfunmat, Ybasis)$fd
  warpfdnames <- unregfd$fdnames
  names(warpfdnames)[3] <- paste("Warped", names(regnames)[1])
  warpfd$fdnames <- warpfdnames
  Ybasis0 <- create.bspline.basis(x0lim, nbasis)
  warpinvfd <- smooth.basis(xfine0, hinvmat, Ybasis0)$fd
  warpfdnames <- unregfd$fdnames
  names(warpfdnames)[3] <- paste("Warped", names(regnames)[1])
  warpinvfd$fdnames <- warpfdnames
  Wfd <- fd(Wcoef, Wbasis)
  return(list(regfd = regfd, warpfd = warpfd, warpinvfd = warpinvfd, 
              Wfd = Wfd))
}

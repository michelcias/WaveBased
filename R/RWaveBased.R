#########################################################################################

PHI <- function(x, J, family = "Daublets", filter.size = 20, prec.wavelet = 30,
                periodic = TRUE){
  
  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry! We don't work with complex data. Only the real part will be
            considered.")
  }
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'.")
  
  PHI <- .Call("C_PHImat", x, as.integer(J), as.integer(fam), as.integer(filter.size),
               as.integer(prec.wavelet), as.integer(periodic))
  return(PHI)
  
}
#########################################################################################

PSI <- function(x, J, family = "Daublets", filter.size = 20, prec.wavelet = 30,
                 periodic = TRUE){
  
  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry! We don't work with complex data. Only the real part will be
            considered.")
  }
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'.")
  
  PSI <- .Call("C_PSImat", x, as.integer(J), as.integer(fam), as.integer(filter.size),
               as.integer(prec.wavelet), as.integer(periodic))
  
  return(PSI)
  
}
#########################################################################################

# wavelet.basis <- function(x, tau = 1, j0 = 0, J, filter.number = 10, family = "DaubLeAsymm",
# prec.wavelet = 30, periodic = TRUE){
#     
#     if(j0 > J) stop("The coarsest level j0 cannot be greater than the finest level J!")
#     
#     if((tau != 1) & (periodic)){
#         tau <- 1
#         warning("The periodic version was calculated with tau = 1.")
#     }
#     
#     if(j0 == J){
#         if(periodic)
#         cnames <- paste("phi[", j0, ",", 0:(tau * 2^j0 - 1), "]", sep = "")
#         
#         else{
#             sup <- support(filter.number, family)
#             phi.lh <- sup$phi.lh
#             phi.rh <- sup$phi.rh
#             kminphi <- c(ceiling(tau * 2^j0 * min(x) - phi.rh + 1e-9))
#             kmaxphi <- c(floor(tau * 2^j0 * max(x) - phi.lh))
#             cnames <- paste("phi[", tau * 2^j0, ",", kminphi:kmaxphi, "]", sep = "")
#         }
#         
#         mat <- PHI(x, tau, j0, filter.number, family, prec.wavelet, periodic)
#         colnames(mat) <- cnames
#         return(mat)
#     }
#     
#     n <- length(x)
#     p <- tau * 2^(j0:(J-1))
#     
#     if(periodic){
#         ncolphi <- p[1]
#         kminphi <- 0
#         kmaxphi <- ncolphi - 1
#         
#         ncolpsi <- p
#         kminpsi <- rep(0, length.out = J-j0)
#         kmaxpsi <- ncolpsi - 1
#     }
#     
#     else{
#         sup <- support(filter.number, family)
#         phi.lh <- sup$phi.lh
#         phi.rh <- sup$phi.rh
#         psi.lh <- sup$psi.lh
#         psi.rh <- sup$psi.rh
#         
#         kminphi <- c(ceiling(p[1] * min(x) - phi.rh + 1e-9))
#         kmaxphi <- c(floor(p[1] * max(x) - phi.lh))
#         ncolphi <- kmaxphi - kminphi + 1
#         
#         kminpsi <- c(ceiling(p * min(x) - psi.rh + 1e-9))
#         kmaxpsi <- c(floor(p * max(x) - psi.lh))
#         ncolpsi <- kmaxpsi - kminpsi + 1
#     }
#     
#     ncolmat <- ncolphi + sum(ncolpsi)
#     
#     mat <- matrix(nrow = n, ncol = ncolmat)
#     cnames <- 1:ncolmat
#     
#     mat[, 1:ncolphi] <- PHI(x, tau, j0, filter.number, family, prec.wavelet, periodic)
#     cnames[1:ncolphi] <- paste("phi[", j0, ",", kminphi:kmaxphi, "]", sep = "")
#     idx <- ncolphi
#     
#     for(j in j0:(J-1)){
#         mat[, idx + 1:ncolpsi[j-j0+1]] <- PSI(x, tau, j, filter.number, family,
#                                               prec.wavelet, periodic)
#         cnames[idx + 1:ncolpsi[j-j0+1]] <- paste("psi[", j, ",", kminpsi[j-j0+1]:kmaxpsi[j-j0+1], "]", sep = "")
#         idx <- idx + ncolpsi[j-j0+1]
#     }
#     
#     colnames(mat) <- cnames
#     return(mat)
# }
#########################################################################################

wbasis <- function(x, j0 = 0, J, family = "Daublets", filter.size = 20, prec.wavelet = 30,
                   periodic = TRUE){
  
  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry! We don't work with complex data. Only the real part will be
            considered.")
  }
  
  if(j0 < 0)
    stop("'j0' must be an integer greater than or equal to zero.")
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'.")
  
  wmat <- .Call("C_WavBasis", x, as.integer(j0), as.integer(J), as.integer(fam), as.integer(filter.size),
               as.integer(prec.wavelet), as.integer(periodic))
  
  return(wmat)
  
}
#########################################################################################

wdensity <- function(data, from, to, length.obs = 250, wf = NULL, power.dens,
                     J0 = 0, J1, family = "Daublets", filter.size = 10, prec.wavelet = 30,
                     dens.biased, warped = TRUE, warp.fun, dwarp.fun, rescale = TRUE,
                     eps = 1.9^(-J1), thresh = "hard", plot = TRUE, main, xlab, ylab,
                     type, ...){
  
  if (!any(is.finite(data))){
    data <- data[is.finite(data)]
    warning("At least one element of 'data' is not finite. Non-finite elements were removed.")
  }
  
  if(eps < 0 | eps > 1)
    stop("'eps' must belong to the interval [0,1].")
  
  thresh <- tolower(substring(thresh, 1, 1))
  if(!thresh %in% c("l","h","s"))
    stop("The options available for 'thresh' are 'linear', 'hard' and 'soft'.")
  
  if(missing(power.dens))
    power.dens <- 0.5
  if(power.dens < 0.5)
    stop("'power.dens' must be greater than or equal to 0.5.")
  
  ry <- range(data)
  sy <- sort(data)
  n <- length(data)
  
  if(missing(from))
    from <- sy[1]
  if(missing(to))
    to <- sy[n]
  if(from < sy[1])
    stop("'from' must be greater than or equal to min(data).")
  if(to > sy[n])
    stop("'to' must be less than or equal to max(data).")
  if(from > to)
    stop("'from' must be smaller than or equal to 'to'.")
  if(missing(dens.biased)){
    dens.bias <- density.default(x = sy, bw = "sj", from = sy[1], to = sy[n])
    dens.biased <- approx(x = dens.bias$x, y = dens.bias$y, xout = sy)$y
  }
  
  y.vals <- seq(from = from, to = to, length.out = length.obs)
  
  if(warped){
    scale <- 1
    if(missing(warp.fun)){
      Hyc <- seq_len(length.out = n)/n
      hy <- dens.biased
      Hy.vals <- approx(x = sy, y = Hyc, xout = y.vals)$y
    }
    else if(!is.function(warp.fun))
      stop("'warp.fun' must be a function.")
    else if(missing(dwarp.fun))
      stop("'dwarp.fun' must be a function.")
    else if(!is.function(dwarp.fun))
      stop("'dwarp.fun' must be a function.")
    else{
      Hyc <- warp.fun(sy)
      hy <- dwarp.fun(sy)
      Hy.vals <- warp.fun(y.vals)
    }
  }
  else if(rescale){
    a <- eps*diff(ry)/(1 - 2*eps)
    location <- ry[1] - a
    scale <- 2*a + diff(ry)
    Hyc <- (sy - location)/scale
    hy <- rep_len(x = 1/scale, length.out = n)
    Hy.vals <- (y.vals - location)/scale
  }
  else{
    Hyc <- sy
    scale <- 1
    hy <- rep_len(x = 1, length.out = n)
    Hy.vals <- y.vals
  }
  
  mbasis <- PHI(x = c(Hyc, Hy.vals), J = J1, family = family, filter.size = filter.size,
                prec.wavelet = prec.wavelet, periodic = TRUE)
  
  coefs <- (scale/mean(1/wf(sy)))^(power.dens)*colMeans(mbasis[seq_len(length.out = n),]*dens.biased^(power.dens-1)*hy/wf(sy)^power.dens)
  
  if(thresh == "h"){
    coefs <- wavedec(x = coefs, j0 = J0, family = family,
                     filter.size = filter.size)
    sig <- mad(coefs[-seq_len(length.out = 2^(J1-1))])
    lambda <- sig*sqrt(2*log(2^(J1-1)))
    d <- coefs[-seq_len(length.out = 2^J0)]
    coefs[-seq_len(length.out = 2^J0)] <- d*(abs(d) > lambda)
    coefs <- waverec(x = coefs, j0 = J0, family = family,
                     filter.size = filter.size)
  }
  if(thresh == "s"){
    coefs <- wavedec(x = coefs, j0 = J0, family = family,
                     filter.size = filter.size)
    sig <- mad(coefs[-seq_len(length.out = 2^(J1-1))])
    lambda <- sig*sqrt(2*log(2^(J1-1)))
    d <- coefs[-seq_len(length.out = 2^J0)]
    coefs[-seq_len(length.out = 2^J0)] <- sign(d)*(abs(d) - lambda)*(abs(d) > lambda)
    coefs <- waverec(x = coefs, j0 = J0, family = family,
                     filter.size = filter.size)
  }
  
  if(!warped & power.dens == 0.5)
    coefs <- coefs/sqrt(sum(coefs^2))
  
  dens.vals <- drop(mbasis[-seq_len(length.out = n),]%*%coefs)^(1/power.dens)/scale
  
  if((1/power.dens) %% 2)
    dens.vals[dens.vals < 0] <- 0
  
  if(plot){
    if(missing(main))
      main <- paste("Density estimates using power =", power.dens)
    if(missing(xlab))
      xlab <- deparse(substitute(data))
    if(missing(ylab))
      ylab <- "Density"
    if(missing(type))
      type <- "l"
    
    plot(x = y.vals, y = dens.vals, type = type, main = main,
         xlab = xlab, ylab = ylab, ...)
  }
  
  return(invisible(list(y = y.vals, dens = dens.vals)))
  
}
#########################################################################################

wavedec <- function(x, j0 = 0, family = "Daublets", filter.size = 20){
  
  if(j0 < 0 | j0 != trunc(j0))
    stop("'j0' must be an integer greater than or equal to zero.")
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets',
         'Symmlets' and 'Coiflets'.")
  
  if(any(!is.finite(x)))
    stop("'x' must be a numeric vector.")
  
  wdec <- .Call("C_WaveDec", as.numeric(x), as.integer(fam), as.integer(filter.size),
                as.integer(j0))
  
  return(wdec)
  
}
#########################################################################################

waverec <- function(x, j0 = 0, family = "Daublets", filter.size = 20){
  
  if(j0 < 0 | j0 != trunc(j0))
    stop("'j0' must be an integer greater than or equal to zero.")
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets',
         'Symmlets' and 'Coiflets'.")
  
  wrec <- .Call("C_WaveRec", as.numeric(x), as.integer(fam), as.integer(filter.size),
                as.integer(j0))
  
  return(wrec)
  
}

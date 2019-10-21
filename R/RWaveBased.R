##########################################################################################

PHI <- function(x, J, family = "Daublets", taps = 20, prec.wavelet = 30,
                periodic = TRUE){
  
  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry! We don't work with complex data. Only the real part will be
            considered.")
  }
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'.")
  
  PHI <- .Call("C_PHImat", x, as.integer(J), as.integer(fam), as.integer(taps),
               as.integer(prec.wavelet), as.integer(periodic))
  return(PHI)
  
}
##########################################################################################

PSI <- function(x, J, family = "Daublets", taps = 20, prec.wavelet = 30,
                 periodic = TRUE){
  
  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry! We don't work with complex data. Only the real part will be
            considered.")
  }
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'.")
  
  PSI <- .Call("C_PSImat", x, as.integer(J), as.integer(fam), as.integer(taps),
               as.integer(prec.wavelet), as.integer(periodic))
  
  return(PSI)
  
}
##########################################################################################

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
##########################################################################################

wbasis <- function(x, j0 = 0, J, family = "Daublets", taps = 20, prec.wavelet = 30,
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
  
  wmat <- .Call("C_WavBasis", x, as.integer(j0), as.integer(J), as.integer(fam), as.integer(taps),
               as.integer(prec.wavelet), as.integer(periodic))
  
  return(wmat)
  
}



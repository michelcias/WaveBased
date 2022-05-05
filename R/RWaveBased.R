################################################################################

PHI <- function(x, J, family = "Daublets", filter.size = 20, prec.wavelet = 30,
                periodic = TRUE, wavelet.filter){
  
  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c", "o"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'. You can also use 'Own', if you provide your wavelet.filter.")
  
  if(missing(wavelet.filter)){
    if(fam == 4)
      stop("Provide your own filter or choose an available family.")
    else
      wavelet.filter <- 0
  }
  else{
    if(any(fam == 1:3))
      warning("wavelet.filter provided, but family is not 'Own'. The arguments family and filter.size were ignored.")
    fam <- 4
  }
  
  PHI <- .Call("C_PHImat", as.double(x), as.integer(J), as.integer(fam),
               as.integer(filter.size), as.integer(prec.wavelet),
               as.integer(periodic), as.double(wavelet.filter))
  return(PHI)
  
}
################################################################################

PSI <- function(x, J, family = "Daublets", filter.size = 20, prec.wavelet = 30,
                 periodic = TRUE, wavelet.filter){
  
  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c", "o"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'. You can also use 'Own', if you provide your wavelet.filter.")
  
  if(missing(wavelet.filter)){
    if(fam == 4)
      stop("Provide your own filter or choose an available family.")
    else
      wavelet.filter <- 0
  }
  else{
    if(any(fam == 1:3))
      warning("wavelet.filter provided, but family is not 'Own'. The arguments family and filter.size were ignored.")
    fam <- 4
  }
  
  PSI <- .Call("C_PSImat", as.double(x), as.integer(J), as.integer(fam),
               as.integer(filter.size), as.integer(prec.wavelet),
               as.integer(periodic), as.double(wavelet.filter))
  
  return(PSI)
  
}
################################################################################

wbasis <- function(x, j0 = 0, J, family = "Daublets", filter.size = 20,
                   prec.wavelet = 30, periodic = TRUE, wavelet.filter){
  
  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }
  
  if(j0 < 0)
    stop("j0 must be a non-negative integer.")
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c", "o"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'. You can also use 'Own', if you provide your wavelet.filter.")
  
  if(missing(wavelet.filter)){
    if(fam == 4)
      stop("Provide your own filter or choose an available family.")
    else
      wavelet.filter <- 0
  }
  else{
    if(any(fam == 1:3))
      warning("wavelet.filter provided, but family is not 'Own'. The arguments family and filter.size were ignored.")
    fam <- 4
  }
  
  wmat <- .Call("C_WavBasis", as.double(x), as.integer(j0), as.integer(J),
                as.integer(fam), as.integer(filter.size),
                as.integer(prec.wavelet), as.integer(periodic),
                as.double(wavelet.filter))
  
  return(wmat)
  
}
################################################################################

wdensity <- function(data, from, to, length.obs = 250, wf = NULL, power.dens,
                     J0 = 0, J1, family = "Daublets", filter.size = 10,
                     prec.wavelet = 30, wavelet.filter, dens.biased, 
                     warped = TRUE, warp.fun, dwarp.fun, rescale = TRUE,
                     eps = 1.9^(-J1), thresh = "hard", plot = TRUE, main, xlab,
                     ylab, type, ...){
  
  if (!any(is.finite(data))){
    data <- data[is.finite(data)]
    warning("At least one element of 'data' is not finite. Non-finite elements were removed.")
  }
  
  if(eps < 0 | eps > 1)
    stop("'eps' must belong to the interval [0,1].")
  
  thresh <- tolower(substring(thresh, 1, 1))
  if(!thresh %in% c("l","h","s"))
    stop("The options available for thresh are 'linear', 'hard' and 'soft'.")
  
  if(missing(power.dens))
    power.dens <- 0.5
  if(power.dens < 0.5)
    stop("power.dens must be greater than or equal to 0.5.")
  
  ry <- range(data)
  sy <- sort(data)
  n <- length(data)
  
  if(missing(from))
    from <- sy[1]
  if(missing(to))
    to <- sy[n]
  
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
      stop("warp.fun must be a function.")
    else if(missing(dwarp.fun))
      stop("dwarp.fun must be a function.")
    else if(!is.function(dwarp.fun))
      stop("dwarp.fun must be a function.")
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
  
  mbasis <- PHI(x = c(Hyc, Hy.vals), J = J1, family = family,
                filter.size = filter.size, prec.wavelet = prec.wavelet,
                periodic = TRUE, wavelet.filter = wavelet.filter)
  
  matvals <- mbasis[seq_len(n),]*dens.biased^(power.dens-1)*hy/wf(sy)^power.dens
  coefs <- (scale/mean(1/wf(sy)))^(power.dens)*colMeans(matvals)
  
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
  
  dens.vals <- drop(mbasis[-seq_len(n),]%*%coefs)^(1/power.dens)/scale

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
################################################################################

wavedec <- function(x, j0 = 0, family = "Daublets", filter.size = 20,
                    wavelet.filter){
  
  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }
  
  if(any(!is.finite(x)))
    stop("'x' must be a numeric vector.")
  
  if(j0 < 0)
    stop("j0 must be a non-negative integer.")
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c", "o"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'. You can also use 'Own', if you provide your wavelet.filter.")
  
  if(missing(wavelet.filter)){
    if(fam == 4)
      stop("Provide your own filter or choose an available family.")
    else
      wavelet.filter <- 0
  }
  else{
    if(any(fam == 1:3))
      warning("wavelet.filter provided, but family is not 'Own'. The arguments family and filter.size were ignored.")
    fam <- 4
  }
  
  wdec <- .Call("C_WaveDec", as.double(x), as.integer(fam),
                as.integer(filter.size), as.integer(j0),
                as.double(wavelet.filter))
  
  return(wdec)
  
}
################################################################################

waverec <- function(x, j0 = 0, family = "Daublets", filter.size = 20,
                    wavelet.filter){
  
  if(is.complex(x)){
    x <- Re(x)
    warning("Sorry, we don't work with complex data. Only the real part was considered.")
  }
  
  if(any(!is.finite(x)))
    stop("'x' must be a numeric vector.")
  
  if(j0 < 0)
    stop("j0 must be a non-negative integer.")
  
  fam <- which(tolower(substring(family, 1, 1)) == c("d", "s", "c", "o"))
  
  if(length(fam) == 0)
    stop("Unknown family. The families available are 'Daublets', 'Symmlets' and 'Coiflets'. You can also use 'Own', if you provide your wavelet.filter.")
  
  if(missing(wavelet.filter)){
    if(fam == 4)
      stop("Provide your own filter or choose an available family.")
    else
      wavelet.filter <- 0
  }
  else{
    if(any(fam == 1:3))
      warning("wavelet.filter provided, but family is not 'Own'. The arguments family and filter.size were ignored.")
    fam <- 4
  }
  
  wrec <- .Call("C_WaveRec", as.double(x), as.integer(fam),
                as.integer(filter.size), as.integer(j0),
                as.double(wavelet.filter))
  
  return(wrec)
  
}

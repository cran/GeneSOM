som <- function(data, xdim, ydim,
                init = "linear",
                alpha = NULL, alphaType = "inverse",
                neigh = "gaussian", topol = "rect",
                radius = NULL, rlen = NULL) {
  SampleInit <- function(xdim, ydim) {
    ind <- sample(1:dim(data)[1], size=xdim*ydim)
    data[ind,]
  }

  RandomInit <- function(xdim, ydim) {
    matrix(rnorm(xdim * ydim * dim(data)[2]), xdim*ydim)
  }

  LinearInit <- function(xdim, ydim) {
    ## get the first two principle components
    pc <- princomp(data)$loadings[,1:2]
    mn <- apply(data, 2, mean)
    ans <- matrix(NA, xdim * ydim, dim(data)[2])
    for (i in 1:(xdim*ydim)) {
      xf <- 4 * (i - 1) %% xdim / (xdim - 1) - 2
      yf <- 4 * (i - 1) %/% xdim / (ydim - 1) - 2
      ans[i, ] <- mn + xf * pc[,1] + yf * pc[,2]
    }
    ans
  }

  INIT <- c("sample", "random", "linear")
  init.type <- pmatch(init, INIT)
  if (is.na(init.type)) stop("Init only supports `sample', `ransom', and `linear'")
  
  ALPHATYPE <- c("linear", "inverse")
  alpha.type <- pmatch(alphaType, ALPHATYPE)
  if (is.na(alpha.type)) stop("AlphaType only supports `linear' and `inverse'.")

  NEIGH <- c("bubble", "gaussian")
  neigh.type <- pmatch(neigh, NEIGH)
  if (is.na(neigh.type)) stop("Neigh only supports `bubble' and `gaussian'")

  TOPOL <- c("rect", "hexa")
  topol.type <- pmatch(topol, TOPOL)
  if (is.na(topol.type)) stop("Topol only supports `rect'.")

  if (init.type == 1) code <- SampleInit(xdim, ydim)
  if (init.type == 2) code <- RandomInit(xdim, ydim)
  if (init.type == 3) code <- LinearInit(xdim, ydim)

  if (is.null(alpha)) alpha <- c(0.05, 0.02)
  else if (length(alpha) == 1) alpha <- c(alpha, alpha/2)
  
  if (is.null(radius)) radius <- c(min(xdim, ydim), min(3, min(xdim, ydim)))
  else if (length(radius) == 1) radius <- c(radius, max(2, radius))

  if (is.null(rlen)) rlen <- c(dim(data)[1] * 2, dim(data)[1] * 10)
  else if (length(rlen) == 1) rlen <- c(rlen, rlen*10)
  
  foo <- .C("som", data=as.double(data),
            as.integer(dim(data)[1]), as.integer(dim(data)[2]),
            code = as.double(code),
            xdim = as.integer(xdim),
            ydim = as.integer(ydim),
            alpha = as.double(alpha),
            alphaType = as.integer(alpha.type),
            neigh = as.integer(neigh.type),
            topol = as.integer(topol.type),
            radius = as.double(radius),
            rlen = as.integer(rlen),
            vis = double(dim(data)[1] * 3))
  visual <- matrix(foo$vis, dim(data)[1], 3)
  visual <- as.data.frame(visual)
  names(visual) <- c("x", "y", "qerror")
  ans <- list(data=data, init=init,
              xdim=xdim, ydim=ydim,
              code=matrix(foo$code, xdim*ydim, dim(data)[2]),
              visual=visual,
              alpha = alpha, alphaType = alphaType,
              neigh = neigh, topol = topol,
              radius = radius, rlen = rlen)
  class(ans) <- "som"
  ans
}

##dyn.load("/home/jyan/work/SOM/GeneSOM/src/test.so")
##foo <- som(matrix(rnorm(100), 25, 4), 3, 4)

filtering <- function(x, lt=20, ut=16000, mmr=3, mmd=200) {
  if (!is.matrix(x)) x <- as.matrix(x)
  ## floor and ceiling
  n <- nrow(x)
  x[x < lt] <- lt
  x[x > ut] <- ut
  tmp1 <- apply(x, 1, max)
  tmp2 <- apply(x, 1, min)
  wch <- (1:n)[ (tmp1/tmp2 > mmr) & (tmp1 - tmp2 > mmd)]
  x[wch,]
}

normalize <- function(x, byrow=TRUE) {
  if (is.vector(x))
    scale(x)
  else if (is.matrix(x) || is.data.frame(x)) {
    if (byrow) t(apply(x, 1, scale))
    else apply(x, 2, scale)
    }
  else stop("The object to be normalized must be vector, matrix, or dataframe.\n")
}

inrange <- function (x, xlim) {
  (x > xlim[1] & x < xlim[2])
}
  
ciplot <- function(x, y, se, n=1, d, ywindow) {
  ind <-  (inrange(y+n*se, ywindow) & inrange(y-n*se, ywindow))
  if (any(ind)) {
    x <- x[ind];y <- y[ind];se <- se[ind]
    segments(x, y+n*se, x, y-n*se)
    segments(x-d, y+n*se, x+d, y+n*se)
    segments(x-d, y-n*se, x+d, y-n*se)
  }
}

plotcell <- function(x, y, dat, code, n, sdbar=1, ylim, yadj) { 
##  yadj <- 0.1
  text(x+1/2, y+(1-yadj/2), paste("n=", n, sep=""))
  if (!is.data.frame(dat)) dat <- as.data.frame(dat)
  ylen <- diff(ylim)
  ## n <- nrow(dat)
  mm <- code
  l <- length(code)
  if (n > 1) {
##      l <- ncol(dat)
##      mm <- sapply(dat, mean)
    ss <- sapply(dat, sd)/sqrt(n)
  }
  else {
##      mm <- dat[,1]
##      l <- length(mm)
##      print(mm)
    ss <- rep(0, l)
  }
  lines(x+((1:l)-1/2)/l, y+1+(ylim[1] + mm)*(1-yadj)/ylen)
  if (sdbar > 0 && n > 1)
    ciplot(x+((1:l)-1/2)/l, y+1+(ylim[1] + mm)*(1-yadj)/ylen,
           ss, n=sdbar, 1/100, ywindow=c(y, y+1-yadj))
}

somgrids <- function(xdim, ydim, color,
                     yadj=0.1, hexa, ntik, ylim) {
  if (color) color <- rainbow(ydim*xdim, start=.7, end=.1)
  else color <- rep("gray", xdim*ydim)
  for (i in 0:(ydim-1)) {
    if (hexa) d <- (i %% 2)/2
    else d <- 0
    lines(c(0+d,xdim+d), c(i,i))
    for (j in 0:xdim) {
      segments(j+d, i, j+d, i+1)
      if (j == xdim) break
      rect(j+d, i+1-yadj, j+1+d, i+1, col=color[j*ydim+i+1])
    }
    lines(c(0+d,xdim+d), c(i+1,i+1))
    if (i %% 2 == 1) axis(2, seq(i, i+1-yadj, length=ntik), seq(ylim[1], ylim[2], length=ntik))
    else axis(4, seq(i, i+1-yadj, length=ntik), seq(ylim[1], ylim[2], length=ntik))
  }
}

plot.som <- function(obj, sdbar=1, ylim=c(-3, 3), color=TRUE, ntik=3, yadj=0.1,
                     xlab="", ylab="", ...) {
 if (class(obj) != "som" ) stop("The funciton must apply to a som object.\n")
 hexa <- (obj$topol == "hexa")
 if (hexa) d <- 1/2
 else d <- 0
 xdim <- obj$xdim; ydim <- obj$ydim
 plot(c(0,xdim+d), c(0,ydim), xlim=c(0,xdim+d), ylim=c(0, ydim),
      type="n", xlab=xlab, ylab=ylab, axes=F, ...)
 axis(1, 0:xdim, 0:xdim)
 
 somgrids(xdim, ydim, color=color, yadj=yadj, hexa=hexa, ylim=ylim, ntik=ntik)
 
 for (i in 0:(ydim-1)) {
   if (hexa) d <- (i %% 2)/2
   else d <- 0
   for (j in 0:(xdim-1)) {
     ind <- obj$visual$x==j & obj$visual$y==i
     n <- length(ind[ind])
     plotcell(j+d, i,
              obj$data[ind, ], obj$code[i*xdim + j+1,], n,
              sdbar=sdbar, ylim=ylim, yadj=yadj)
   }
 }
}

print.som <- function(obj) {
  tmp <- summary(obj)
  print(tmp)
}

summary.som <- function(obj) {
  tmp <- list(init=obj$init, neigh=obj$neigh, topol=obj$topol,
              xdim=obj$xdim, ydim=obj$ydim,
              alpha=obj$alpha, alphaType=obj$alphaType,
              radius=obj$radius, rlen=obj$rlen)
  tmp
}

summary.print.som <- function(obj) {
  print(summary(tmp))
}

qerror <- function(obj, radius=1) {
  if (class(obj) != "som" ) stop("The funciton must apply to a som object.\n")
  ALPHATYPE <- c("linear", "inverse")
  alpha.type <- pmatch(obj$alphaType, ALPHATYPE)
  if (is.na(alpha.type)) stop("AlphaType only supports `linear' and `inverse'.")

  NEIGH <- c("bubble", "gaussian")
  neigh.type <- pmatch(obj$neigh, NEIGH)
  if (is.na(neigh.type)) stop("Neigh only supports `bubble' and `gaussian'")

  TOPOL <- c("rect", "hexa")
  topol.type <- pmatch(obj$topol, TOPOL)
  if (is.na(topol.type)) stop("Topol only supports `rect'.")
  ##data <- obj$data
  ##code <- obj$code
  ##visual <- as.matrix(obj$visual)
  ans <- .C("find_qerror", as.double(obj$data),
     as.integer(dim(obj$data)[1]), as.integer(dim(obj$data)[2]),
     as.double(obj$code),
     as.integer(obj$xdim), as.integer(obj$ydim),
     as.integer(alpha.type),
     as.integer(neigh.type), as.integer(topol.type),
     as.double(radius), 
     as.double(as.matrix(obj$visual)), qerror=double(1))$qerror
  ans / dim(obj$data)[1]
}

# MAIN FUNCTIONS FOR SPECTRAL CHANGE POINT ESTIMATION #####


## PART 1: SPECTRAL BLOCK ESTIMATION AND PROJECTION ####

#' Bartlett weights.
#' @keywords internal
Bartlett.weights <- function(x) {
  1 - abs(x)
}

#' Function to apply the normal density transformation.
#' This function is helpful to relieve the outlier/heavy tail effect.
#' @param x a vector
#' @returns a transformed vector
#' @keywords internal
NormalTrans <- function(x) {
  # x is a vector
  len <- length(x)
  fn <- ecdf(x)
  xtrans <- qnorm(fn(x) - 1 / 2 / len)
  return(xtrans)
}


#' Spectral matrix estimate in blocks
#'
#' Lag window estimate of spectral matrix in non-overlapping blocks.
#'
#' @param y A `T` by `m` matrix for a multivariate time series of dimension `m` and length `T` .
#' @param L The block size, a scalar `0<L<T`. Could be set subject to user preference. At least 50 for a good estimate of spectrum within block.
#' @param freq The frequencies on which the spectral density. will be estimated, can be specified due to expert knowledge or user interest. The default values are 18 values evenly distributed between 0 and pi.
#' @param R (optional) The bandwidth, a scalar `0<R<L` with default `R= floor(L^{1/3})`.
#' @param band (optional) Use it when users want to assess frequency bands instead of individual frequencies. An integer vector indicating the number in each band. Specifically, the vector is of length \eqn{a} and sum to \eqn{b}, where \eqn{a} is the number of groups, and \eqn{b}  equals to length of freq. Default value means no band.
#' @param scale (optional, default `TRUE`) Logical indicating whether marginal scale will be implemented before estimation.
#' @param normaltrans (optional, default `TRUE`) Logical indicating whether a marginal nonparametric normal transformation will be implemented before estimation.
#' @param coherency (optional, default `FALSE`)  Use it if users want to assess the coherency matrix, instead of spectral matrix. If `TRUE`, the spectral matrix will be normalized such that it has unit diagonal.
#' @returns Returns the spectral estimate, which is a four-dimensional array, of dim by dim by freq (or band) by block
#' @keywords internal
SpecBlock <- function(y,
                      L,
                      freq = seq(0, pi, length.out=19)[-1],
                      R = floor(L ^ (1 / 3)),
                      band = rep(1, length(freq)),
                      scale = TRUE,
                      normaltrans = TRUE,
                      coherency = FALSE) {
  if (is.vector(y))
    y <- matrix(y, ncol = 1)
  T <- dim(y)[1]
  m <- dim(y)[2]
  if (normaltrans)
    y <- apply(y, 2, NormalTrans)
  if (scale) {
    ys <- scale(y)
  } else {
    ys <- y
  }

  B <- floor(T / L)  # only keep complete blocks

  w <- Bartlett.weights(((-R):R) / R)

  nfreq <- length(freq)
  nband <- length(band)
  band.ind <- rep(1:nband, times = band)

  bs <- bs_block(ys = ys, w = w, freq = freq, band = band, band_ind = band.ind, R = R, B = B, L = L, coherency = coherency)
  bs <- array(bs, dim = c(m, m, nband, B))
  bs <- Re(bs)

  if (!is.null(dimnames(y)[[2]])) {
    dimnames(bs)[[1]] <- dimnames(bs)[[2]] <- dimnames(y)[[2]]
  } else {
    dimnames(bs)[[1]] <- dimnames(bs)[[2]] <- paste("V", 1:m, sep = "")
  }
  dimnames(bs)[[3]] <- round(tapply(freq, band.ind, mean), 3)
  dimnames(bs)[[4]] <- paste("B", 1:B, sep = "")
  bs
}

#' Piecesewise block mean of spectral matrix estimate
#'
#' This function is used to show the piecewise mean after estiamting the change points.
#' @inheritParams SpecBlock
#' @param bs The block spectral, if not given, will be estimated from data using `SpecBlock`.
#' @param cp The estimated cp location in block scale.
#' @returns Returns the piecesewise block mean of spectral estimate, which is a four-dimensional array, of dim by dim by freq (or band) by number of stationary periods.
#' @keywords internal
SpecBlockMean <- function(y = NULL,
                          L = NULL,
                          freq = seq(0, pi, length.out=19)[-1],
                          R = floor(L ^ (1 / 3)),
                          band = rep(1, length(freq)),
                          scale = TRUE,
                          normaltrans = TRUE,
                          coherency = FALSE,
                          bs = NULL,
                          cp = NULL){
  if(missing(bs)){
    bs <- SpecBlock(
      y = y,
      L = L,
      R = R,
      freq = freq,
      band = band,
      scale = scale,
      normaltrans = normaltrans,
      coherency = coherency
    )
  }

  B <- dim(bs)[[4]]
  len <- diff(c(0, cp, B))
  ind <- rep(1:length(len), len)

  mean <- apply(bs, c(1, 2, 3), function(x) {
    tapply(x, ind, mean)
  })
  mean <- aperm(mean, c(2, 3, 4, 1))
  return(mean)

}


#' Function for the spectral tensor CUSUM, designed specifically for the symmetric spectral matrix.
#' @param x 4-order tensor of dim by dim by freq (or band) by block (at)
#' @returns 4-order CUSUM tensor of dim by dim by freq (or band) by block (at)-1
#' @keywords internal
TensorCusum <- function(x) {
  m <- dim(x)[1]
  nfreq <- dim(x)[3]
  bc0 <- array(0, dim = c(dim(x)[1:3], dim(x)[4] - 1))
  for(k in 1:nfreq){
    if (m == 1){
      bc0[, , k, ] <- - cusum_vec(x[, , k, ])
    } else {
      bc0[, , k, ] <- - cusum_symm_tensor(x[, , k, ])
    }
  }
  return(bc0)
}

#' Function for frequency-specific projection
#' @param bs.part the four order spectral tensor, of dim by dim by freq (or band) by block (or at)
#' @param sparsity the sparsity value
#' @param initial the initialization method, can be either aggregation or unfolding.
#' @returns a list including information of the projection
#' @keywords internal
SpecProj <- function(bs.part, sparsity = NULL, initial = "aggregation") {
  m <- dim(bs.part)[1]
  len <- dim(bs.part)[4]
  bc0 <- TensorCusum(bs.part)

  nfreq <- dim(bs.part)[3]

  bcp.l <- matrix(0, nfreq, len - 1)
  rownames(bcp.l) <- dimnames(bs.part)[[3]]
  pro.l <- matrix(0, nfreq, m)
  rownames(pro.l) <- dimnames(bs.part)[[3]]
  colnames(pro.l) <- dimnames(bs.part)[[1]]
  for (i in 1:nfreq) {
    if (m > 1) {
      # multivariate, projection needed
      result <- mystensor(bc0[, , i, ], sparsity = sparsity, initial = initial)
      pro.l[i, ] <- result$a
      pro <- matrix(result$a, nrow = 1)
      # projected spectral density
      bsp <- as.vector(Re(
        ttl_rcpp(tnsr = bs.part[, , i, ], list_mat = list(pro, pro), ms = c(1, 2))))
      # projected cusum
      bc0p <- Re(cusum_vec(bsp))
      # projected cusum normalized by mean of spectral density
      bcp.l[i, ] <- abs(bc0p / mean(bsp))

    } else {
      # univariate, projection not needed
      bcp.l[i, ] <- abs(bc0[, , i, ] / mean(bs.part[, , i, ]))
    }

  }
  cpseries.l <- apply(pro.l, 1, function(x) {
    which(x != 0)
  })
  if (!is.matrix(cpseries.l))
    cpseries.l <- matrix(cpseries.l, nrow = 1)
  cpseries.l <- t(cpseries.l)
  if (length(cpseries.l > 0)) {
    rownames(cpseries.l) <- dimnames(bs.part)[[3]]
    colnames(cpseries.l) <- NULL
  }


  output <- NULL
  output$bcp.l <- bcp.l
  output$cpseries.l <- cpseries.l
  output$pro.l <- pro.l
  output$sparsity <- sparsity
  return(output)
}




## PART 2: THRESHOLD COMPUTATION ####

#' Function for compute the list of maximum projected cusum  values based on many frequency bootstrap samples, in order to get the threshold
#' @param bs the block spectral estimate by function `SpecBlock`, of dim by dim by freq (or band) by block
#' @param sparsity  the sparsity value
#' @param Nb the number of bootstraps for obtaining the threshold value
#' @param do.parallel number of copies of R running in parallel
#' @param quan take the \eqn{quan}-th quantile among bootstraps for the threshold value
#' @param initial the initialization method, can be either aggregation or unfolding.
#' @returns a list containing the threshold in all bootstraps, and the final threshold
#' @keywords internal
Tau.compu <- function(bs = NULL,
                      sparsity,
                      Nb = 100,
                      do.parallel = 0,
                      quan = 0.975,
                      initial = "aggregation") {
  m <- dim(bs)[[1]]
  B <- dim(bs)[[4]]

  if (do.parallel > 0) {
    # parallel computation
    future::plan(future::multisession, workers = do.parallel)
    `%mydo%` <- doFuture::`%dofuture%`
  } else {
    `%mydo%` <- foreach::`%do%`
  }


  block.ind <- 1:B
  # robust step, remove blocks with 10% high spec norm in multivariate case
  if (m > 1) {
    block.sum <- apply(apply(bs, c(3, 4), norm, "2"), 2, sum)
    block.rem <- order(block.sum, decreasing = T)[1:(floor(B * 0.1))]
    block.ind <- (1:B)[-block.rem]
  }



  Tau <- function(x) {
    boot.ind <- sample(block.ind, size = B, replace = TRUE)
    bs.boot <- bs[, , , boot.ind, drop = FALSE]
    pro <- SpecProj(bs.boot, sparsity, initial = initial)
    if (m > 1) {
      # in multivariate case, frequency specific threshold
      res <- apply(pro$bcp.l, 1, max)
    } else {
      # in univariate case, overall threshold
      res <- max(pro$bcp.l)
    }
    return(res)
  }


  res <- foreach::foreach(
    x = 1:Nb,
    .combine = "c",
    .options.future = list(
      seed = TRUE,
      packages = c(
        "Rcpp",
        "RcppArmadillo",
        "SpecCp",
        "foreach",
        "future",
        "doFuture"
      )
    )
  ) %mydo% {
    Tau(x)
  }
  res <- matrix(res, ncol = Nb)

  tau <- apply(res, 1, quantile, probs = quan)
  return(list(tau.l = res, tau = tau))
}

## PART 3: MULTIPLE CHANGE POINT DETECTION BY WSBS ####

#' Function to find the maximum point of a positive vector, considering only points around which an interval of length  \eqn{2\delta} are all positive
#' @param stat the positive vector
#' @param delta let the picked location within a interval with positive projected cusum, of length \eqn{2\delta}
#' @returns a list containing the maximum point and the maximum value
#' @keywords internal
find.b <- function(stat, delta) {
  FLAG <- FALSE
  b <- 1
  test.stat <- 0
  stat.update <- stat

  while (!FLAG && sum(stat.update) > 0) {
    b <- which.max(stat.update)
    int <- max(b - delta, 1):min(b + delta, length(stat))
    if (sum(stat[int] == 0) > 0) {
      # if the interval are not all positive, let the point be zero, and find the next maximum one
      stat.update[b] <- 0
    } else {FLAG <- TRUE}
    test.stat <- max(stat.update)
  }
  return(list(b = b, test.stat = test.stat))
}
#'  Function of sparsified sum over frequencies for multivariate time series
#'  @param output the output by `SpecProj`, if not given, will be computed with `bs.part`
#'  @param bs.part the four order spectral tensor, of dim by dim by freq (or band) by block
#'  @param m the dimension, since univariate and multivariate cases are different that the former does not need projection
#' @param sparsity the sparsity value
#' @param tau the threshold
#' @param trim the boundary removal \eqn{\nu}
#' @param initial the initialization method, can be either aggregation or unfolding.
#' @returns a list including the maximum cusum, the maximum location, and informations of the series and frequency.
#' @keywords internal
sparsum <- function(output = NULL,
                    bs.part = NULL,
                    m,
                    sparsity = NULL,
                    tau,
                    trim = NULL,
                    initial = "aggregation") {
  if (missing(output))
    output <- SpecProj(bs.part, sparsity, initial)
  len <- dim(output$bcp.l)[2] + 1
  bcp.l <- output$bcp.l
  if (m > 1) {
    cpseries.l <- output$cpseries.l
    pro.l <- output$pro.l
    tau.m <- matrix(rep(tau, len - 1), ncol = len - 1) # frequency specific threshold
    over <- bcp.l > tau.m
    over.sum <- apply(over, 1, sum)
    cpseries.l[which(over.sum == 0), ] <- 0
  } else if (m == 1) {
    over <- bcp.l > tau # universal threshold
  }


  bcsp <- colSums(bcp.l * over)  #block cusum summary
  bcsp[c(1:(1 + trim), (len - 1 - trim):(len - 1))] <- 0


  search.res <- find.b(stat = abs(bcsp), delta = max(1, floor(trim / 4)))
  cusum <- search.res$test.stat
  b <- search.res$b
  output <- NULL
  output$cusum <- cusum
  output$b <- b
  cpfreq <- which(over[, b] == TRUE)
  if (m > 1) {
    cpseries <- cpseries.l[cpfreq, ]
    pro <- pro.l[cpfreq, ]
    cpfs <- list(
      cpfreq = cpfreq,
      #cpseries = cpseries,
      #pro = pro,
      pro.allfreq = pro.l
    )
  } else {
    cpfs <- list(cpfreq = cpfreq)
  }

  output$cpfs <- cpfs
  return(output)
}

#' Function of wild sparsified binary segmentation
#' @param bs the four order spectral tensor, of dim by dim by freq (or band) by block
#' @param sparsity the sparsity value
#' @param tau the threshold
#' @param trim the boundary removal \eqn{\nu}
#' @param M the number of random intervals in the wild binary segmentation
#' @param do.parallel number of copies of R running in parallel
#' @param single logical, if TRUE, only a single change point will be estimated, and binary segmentation will not be implemented.
#' @param initial the initialization method, can be either aggregation or unfolding.
#' @returns a list including the change point location information, and the change point frequency and series information
#' @keywords internal
wsbs <- function(bs,
                 sparsity,
                 tau,
                 trim = NULL,
                 M = NULL,
                 do.parallel = 0,
                 single = FALSE,
                 initial = " aggregation") {
  B <- dim(bs)[4]
  m <- dim(bs)[1]
  if (missing(M))
    M <- 2 * B
  rnd1 <- sample(1:(B - 2 * trim), M, replace = TRUE)
  rnd2 <- sample(1:(B - 2 * trim), M, replace = TRUE)
  window_s <- pmin(rnd1, rnd2)
  window_e <- pmax(rnd1, rnd2) + 2 * trim

  Getcp_int <- function(x) {
    s_x <- window_s[x]
    e_x <- window_e[x]
    obj <- SpecProj(bs[, , , (s_x):(e_x), drop = FALSE], sparsity, initial)
    return(obj)
  }

  if(M > 0){
    if (do.parallel > 0) {
      future::plan(future::multisession, workers = do.parallel)
      `%mydo%` <- doFuture::`%dofuture%`
    } else {
      `%mydo%` <- foreach::`%do%`
    }
    int.l <- foreach::foreach(x = 1:M,
                              .options.future = list(
                                seed = TRUE,
                                packages = c(
                                  "Rcpp",
                                  "RcppArmadillo",
                                  "SpecCp",
                                  "foreach",
                                  "future",
                                  "doFuture"
                                )
                              )) %mydo% {
                                Getcp_int(x)
                              }

    max.l <- rep(0, M)
    cp.l <- rep(0, M)
    cpfs.l <- NULL
    for (x in 1:M) {
      res <- sparsum(
        output = int.l[[x]],
        m = m,
        sparsity = sparsity,
        tau = tau,
        trim = trim,
        initial = initial
      )

      max.l[x] <- res$cusum
      cp.l[x] <- res$b + window_s[x] - 1
      cpfs.l[[x]] <- res$cpfs
    }
  }


  # View(cbind(window_s,window_e,max.l,cp.l))

  BinSeg <- function(bs, s, e, depth, single) {
    if (e - s <= 2 * trim)
      return(NULL)
    ind <- which(window_s >= s & window_e <= e)
    max.val <- -1
    cp <- 0
    s.int <- s
    e.int <- e

    # the interval with s0 and e0
    obj <- sparsum(
      bs.part = bs[, , , s:e, drop = FALSE],
      m = m,
      sparsity = sparsity,
      tau = tau,
      trim = trim,
      initial = initial
    )


    max.val <- obj$cusum
    cp <- s - 1 + obj$b
    cpfs <- obj$cpfs
    # compare the random intervals with [s0,e0]
    if (sum(ind) > 0) {
      if (max.val < max(max.l[ind])) {
        max.val <- max(max.l[ind])
        ind.max <- ind[which.max(max.l[ind])]
        cp <- cp.l[ind.max]
        cpfs <- cpfs.l[[ind.max]]
        s.int <- window_s[ind.max]
        e.int <- window_e[ind.max]

      }
    }
    ret <- NULL
    ret$loc.block <- cp
    ret$max.proj.cusum <- max.val
    ret$depth <- depth
    ret$s.int <- s.int
    ret$e.int <- e.int
    ret$cpfs <- cpfs
    if (ret$max.proj.cusum <= 0) {
      return(NULL)
    } else if (!single) {
      return(cbind(
        BinSeg(bs, s, cp, depth + 1, single = FALSE),
        ret,
        BinSeg(bs, cp + 1, e, depth + 1, single = FALSE)
      ))
    } else if (single) {
      return(cbind(ret))
    }
  }
  ret <- NULL
  changepoints <- BinSeg(bs, 1, B, depth = 1, single = single)
  cp <- t(matrix(as.numeric(changepoints[-(6), ]), nrow = 5))
  colnames(cp) <- c("loc.block", "max.proj.cusum", "depth", "s.int", "e.int")
  cpfs <- changepoints[6, ]
  if (!is.null(cpfs))
    names(cpfs) <- 1:length(cpfs)
  ret$cp <- cp
  ret$cpfs <- cpfs
  return(ret)
}

## PART 4: CHANGE POINT LOCATION REFINEEMENT ####


#' Function of lag window estimate for projected spectral at a consecutive period this is for the refinement step
#' @param y  the time series, multivariate (\eqn{T*m})
#' @param L the block size, a scalar \eqn{0<L<T}
#' @param R the bandwidth, a scalar \eqn{0<R<L} with default \eqn{R= floor(L^(1/3))}
#' @param freq the frequencies on which the spectral density will be estimated, can be specified due to expert knowledge or user interest
#' @param band an integer vector when frequency band average is implemented, indicating the number in each band. Specifically, the vector is of length \eqn{a} and sum to \eqn{b}, where \eqn{a} is the number of groups, and \eqn{b}  equals to length of freq. Default value means no band.
#' @param scale logical indicating whether marginal scale will be implemented before estimation
#' @param normaltrans logical indicating whether a marginal nonparametric normal transformation will be implemented before estimation
#' @param coherency logical indicating whether a coherency matrix, instead of a spectral matrix will be estimated. If TRUE,we scale the spectral matrix such that it has unit diagonal
#' @param at specifying the time points at which the spectral will be estimate, should be consecutive points, and the minimum is larger than \eqn{L}
#' @param projections  a matrix of dim length(band) by m, which will project the spectral to a of 1 by 1 by freq (or band) by at tensor.

#' @returns returns the projected spectral estimate at `at`, which is a four-dimensional array, of 1 by 1 by freq (or band) by at
#' @keywords internal
SpecBlockOverlapProject <- function(y,
                                    L,
                                    R = floor(L ^ (1 / 3)),
                                    freq,
                                    band = rep(1, length(freq)),
                                    scale = TRUE,
                                    normaltrans = FALSE,
                                    coherency = FALSE,
                                    at = NULL,
                                    projections = NULL) {
  if (is.vector(y))
    y <- matrix(y, ncol = 1)
  T <- dim(y)[1]
  m <- dim(y)[2]
  if (normaltrans)
    y <- apply(y, 2, NormalTrans)
  if (scale) {
    ys <- scale(y)
  } else {
    ys <- y
  }


  B <- floor(T / L)  # only keep complete blocks

  w <- Bartlett.weights(((-R):R) / R)

  nfreq <- length(freq)
  nband <- length(band)
  band.ind <- rep(1:nband, times = band)

  bs <- bs_consecutive(ys = ys, w = w, freq = freq, band = band, band_ind = band.ind, at = at, projections = projections, R = R, L = L, coherency = coherency)
  bs <- array(bs, dim = c(1, 1, nband, length(at)))
  bs <- Re(bs)

  dimnames(bs)[[1]] <- dimnames(bs)[[2]] <- "projected"
  dimnames(bs)[[3]] <- tapply(freq, band.ind, mean)
  dimnames(bs)[[4]] <- at
  bs
}

#' Function to get a refined change point location using overlap blocks, after got change point block locations. The parameters include some parameters in function `SpecBlockOverlapProject`, as well as
#' @param res an internal result by the function `wsbs`
#' @param bs the non-overlapped block estimate
#' @param tau the threshold
#' @param trim the boundary removal
#' @returns the refined change point location
#' @keywords internal
SpecCpOverlap <- function(y,
                          L,
                          R = floor(L ^ (1 / 3)),
                          freq,
                          band = rep(1, length(freq)),
                          scale = TRUE,
                          normaltrans = FALSE,
                          coherency = FALSE,
                          res,
                          bs,
                          tau,
                          trim) {
  if (is.vector(y))
    y <- matrix(y, ncol = 1)
  T <- dim(y)[1]
  m <- dim(y)[2]
  nband <- length(band)
  cp.overlap <- rep(0, nrow(res$cp))

  for (i in 1:nrow(res$cp)) {
    cp.block <- res$cp[i, "loc.block"]
    if (m == 1) {
      cp.pro <- matrix(1, nband, 1)
    } else {
      cp.pro <- res$cpfs[[i]]$pro.allfreq
    }

    # overlapped in the two blocks
    at <- ((L * (cp.block - 1) + 1):(L * (cp.block + 5)))
    aroundbsp <- SpecBlockOverlapProject(
      y = y,
      L = L * 4,
      R = R,
      freq = freq,
      band = band,
      scale = scale,
      normaltrans = normaltrans,
      coherency = coherency,
      at = at,
      projections = cp.pro
    )
    di <- apply(aroundbsp, 3, diff, lag = L * 4)
    di.max <- apply(abs(di), 1, max)
    cp.overlap[i] <- which.max(di.max) + at[1] - 1
  }
  return(cp.overlap)
}




## PART 5: MAIN FUNCTION FOR MULTIVARIATE FREQUENCY CHANGE POINT ESTIMATION ####

#' Change point estimation for multivariate time series from frequency domain
#'
#' Estimate the change point locations using the Spectral method proposed in Zhang and Chan (2025). Also outputs the change point frequency and series information.
#'
#' @param y A `T` by `m` matrix for a multivariate time series of dimension `m` and length `T` .
#' @param L The block size, a scalar `0<L<T`. Could be set subject to user preference. At least 50 for a good estimate of spectrum within block.
#' @param freq The frequencies on which the spectral density. will be estimated, can be specified due to expert knowledge or user interest. The default values are 18 values evenly distributed between 0 and pi.
#' @param R (optional) The bandwidth, a scalar `0<R<L` with default `R= floor(L^{1/3})`.
#' @param band (optional) Use it when users want to assess frequency bands instead of individual frequencies. An integer vector indicating the number in each band. Specifically, the vector is of length \eqn{a} and sum to \eqn{b}, where \eqn{a} is the number of groups, and \eqn{b}  equals to length of freq. Default value means no band.
#' @param sparsity (optional) The sparsity value. Default values are estimated from the data based on the data using one-by-one approach, so users do not need to specify unless they want to override the defaults.
#' @param refine (optional, default `TRUE`) Logical indicating whether a refinement step will be indicated by scanning the localized blocks around the change point detected for multivariate time series
#' @param do.parallel (optional, default 0) Number of copies of R running in parallel
#' @param initial the initialization method, can be either aggregation or unfolding, but aggregation is recommended since it is faster and more stable.
#' @param mean.output (optional, default `FALSE`) Logical indicating whether the piecewise mean spectral will be output
#' @param tau (optional) The threshold. Default values are estimated from the data, so users do not need to specify unless they want to override the defaults.
#' @param M (optional) The number of random intervals in the wild binary segmentation.
#' @param Nb (optional) The number of bootstraps for obtaining the threshold value.
#' @param scale (optional, default `TRUE`) Logical indicating whether marginal scale will be implemented before estimation.
#' @param normaltrans (optional, default `TRUE`) Logical indicating whether a marginal nonparametric normal transformation will be implemented before estimation.
#' @param coherency (optional, default `FALSE`)  Use it if users want to assess the coherency matrix, instead of spectral matrix. If `TRUE`, the spectral matrix will be normalized such that it has unit diagonal.
#' @param single (optional, default `FALSE`) If `TRUE`, only a single change point will be estimated, and binary segmentation will not be implemented. Should be omitted and only be used for a fast conclusion of breaks existence.
#' @param message (optional, default `TRUE`) If `TRUE`, print out the computation process.
#' @returns A list with the following elements
#' \describe{
#'  \item{cp}{A matrix with `k` rows, where `k` is the number of detected change points. Each row include: `loc.block` is the change point block index , `max.proj.cusum` is the maximum projected cusum value at this change point, `depth` is the depth in the wild binary segmentation, `s.int` and `e.int` are the start and end of the interval on which this change point is detected, `loc` is the change point time index `loc.block * L`, `loc.refine` is the refined change point.}
#'  \item{cpfs}{A list with `k` sub-lists, one for each change point, in the order of `cp`. Each sub-list include two elements: `cpfreq` is the change point frequency at this change point, and `pro.allfreq` is a matrix with `length(freq)` (or `length(band)` subject `band` argument setting) rows and `m` columns representing the projection at this change point.}
#'  \item{sparsity}{The specified or data-driven sparsity.}
#'  \item{mean}{The piecewise mean of the spectral estimate (or coherency subject to `coherency` argument setting), a 4-dimensional order of `m` by `m` by `length(freq)` by `k+1`. Could be scaled subject to argument setting.}
#' }
#' @examples
#' \dontrun{
#' # Data generation and change point detection
#' set.seed(1)
#' y <-DGP_mcp(TT = 6000, m = 80,
#' sparsity = 8, q = 4, type = "vma")
#' res <- SpecCp(y = y, L = 75)
#'
#' # Results
#' res$cp # cp location
#' res$cpfs$`1`$cpfreq # cp frequency
#' mat <- res$cpfs$`1`$pro.allfreq
#' CpFreq <- rep(0, nrow(mat))
#' CpFreq[res$cpfs$`1`$cpfreq] <- 1
#' annotation_row <- data.frame(CpFreq = CpFreq)
#' rownames(annotation_row) <- rownames(mat)
#' library(pheatmap)
#' pheatmap(mat, annotation_row = annotation_row,
#' main = "Projection", cluster_rows = FALSE,
#' cluster_cols = FALSE,) # cp frequency-series
#' }
#' @references Zhang and Chan. (2025). Spectral change point estimation for high dimensional time series by sparse tensor decomposition.
#' @author Xinyu Zhang
#' @export
SpecCp <- function(y,
                   L,
                   freq = seq(0, pi, length.out=19)[-1],
                   R = max(1, floor(L ^ (1 / 3))),
                   band = rep(1, length(freq)),
                   sparsity = NULL,
                   refine = TRUE,
                   do.parallel = 0,
                   initial = "aggregation",
                   mean.output = FALSE,
                   tau = NULL,
                   M = NULL,
                   Nb = 100,
                   scale = TRUE,
                   normaltrans = TRUE,
                   coherency = FALSE,
                   single = FALSE,
                   message = TRUE) {
  if (is.vector(y))
    y <- matrix(y, ncol = 1)
  T <- dim(y)[1]
  m <- dim(y)[2]
  B <- floor(T / L)
  trim <- max(1, floor((log(T * m) * B) ^ (2 / 3) / 15))
  if (missing(M))
    M <- 2 * B

  # sparsity
  if (!missing(sparsity)) {
    if (sparsity <= 0 | sparsity != round(sparsity) |
        sparsity > m) {
      warning("Sparsity is not a reasonable integer.")
      return(NULL)
    }
  } else if (m > 1) {
    if(message) print("Sparsity estimation: started...")
    res.one <- list(NULL)
    for (j in 1:m) {
      res <- SpecCp(
        y = y[, j],
        L = L,
        R = R,
        freq = freq,
        band = band,
        refine = FALSE,
        do.parallel = do.parallel,
        initial = initial,
        M = M,
        scale = scale,
        normaltrans = normaltrans,
        coherency = coherency,
        single = TRUE,
        message = FALSE
      )
      res.one[[j]] <- as.numeric(res$cp[, 1])
    }
    sparsity <- sum(lengths(res.one) != 0)
    if(message) print("Sparsity estimation: completed.")
  }


  if(message) print("Spectral matrix estimation: started...")
  bs <- SpecBlock(
    y = y,
    L = L,
    R = R,
    freq = freq,
    band = band,
    scale = scale,
    normaltrans = normaltrans,
    coherency = coherency
  )
  if(message) print('Spectral matrix estimation: completed.')

  # threshold
  if (missing(tau)) {
    if(message) print("Threshold computation: started...")
    quan <- 0.975
    tau <- Tau.compu(
      bs = bs,
      sparsity = sparsity,
      Nb = Nb,
      do.parallel = do.parallel,
      quan = quan,
      initial = initial
    )$tau
    if(message) print('Threshold computation: completed.')
  }

  if(message) print("Block change point location estimation: started...")
  res <- wsbs(
    bs = bs,
    sparsity = sparsity,
    tau = tau,
    trim = trim,
    M = M,
    do.parallel = do.parallel,
    single = single,
    initial = initial
  )
  res$cp <- cbind(res$cp, loc = res$cp[,1] * L)
  res$sparsity <- sparsity
  if(message) print("Block change point location estimation: completed.")

  # piece wise mean
  if (mean.output & nrow(res$cp) > 0){
    mean <- SpecBlockMean(bs = bs, cp = res$cp[, 1])
    res$mean <- mean
  }



  # overlap refinement
  if (refine & nrow(res$cp) > 0 & m>1) {
    if(message) print("Change point location refinement: started...")
    cp.overlap <- SpecCpOverlap(
      y = y,
      L = L,
      R = R,
      freq = freq,
      band = band,
      scale = scale,
      normaltrans = normaltrans,
      coherency = coherency,
      res = res,
      bs = bs,
      tau = tau,
      trim = trim
    )
    res$cp <- cbind(res$cp, loc.refine = cp.overlap)
    if(message) print("Change point location refinement: completed.")
  }


  return(res)
}


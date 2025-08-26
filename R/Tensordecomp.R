
# FUNCTIONS FOR SPARSE TENSOR DECOMPOSTION  #


#' Function of getting the sparse leading eigen vector
#' @param X symmetric matrix
#' @param sparsity sparsity value
#' @param v1.init the initial vector
#' @keywords internal
mytpower <- function(X,
                    sparsity,
                    v1.init,
                    niter = 10,
                    lambda.diff.threshold = 1e-06,
                    trace = FALSE) {
  if (sum(abs(X)) == 0) {
    return(rep(0, dim(X)[1]))
  } else {
    p <- nrow(X)
    if (missing(v1.init)) {
      v1.init <- powerIteration0(
        X = X,
        k = p,
        max_iter = 10,
        lambda_diff_threshold = lambda.diff.threshold,
        trace = trace
      )$v1
    }
    power.iter.out <- powerIteration0(
      X = X,
      k = sparsity,
      v1_init = v1.init,
      max_iter = niter,
      lambda_diff_threshold = lambda.diff.threshold,
      trace = trace
    )
    return(as.vector(power.iter.out$v1))
  }
}


#'Function to get sparse leading term of a partial symmetric (at mode 1 and 2) partial constant (at mode 3) tensor via power method
#' @param T the three order tensor
#' @param a0 the initial vector at mode 1
#' @param sparsity sparsity value at mode 1
#' @param niter maximum number of iterations
#' @keywords internal
myalternative <- function(T, a0, sparsity = NULL, niter = 10) {
    d1 <- dim(T)[1]
    d3 <- dim(T)[3]
    if (is.null(sparsity)) { # if not specified, use that of the initial vector
        sparsity <- sum(a0 != 0)
    }

    a.old <- rep(0, d1)
    c.old <- rep(0, d3)
    a.new <- a0
    c.new <- c.old
    iter.t <- 0
    tt <- T

    EPS <- 1e-04
    while (max(
      min(norm(a.old - a.new, type = "2"), norm(a.old + a.new, type = "2")),
      min(norm(c.old - c.new, type = "2"), norm(c.old + c.new, type = "2"))) >= EPS
           && iter.t < niter) {
        iter.t <- iter.t + 1
        a.old <- matrix(a.new, nrow = 1)
        c.old <- matrix(c.new, nrow = 1)

        c.tmp <- ttl_rcpp(tnsr = tt, list_mat = list(a.old, a.old), ms = c(1, 2))
        if (sum(c.tmp) < 0) # unify the signs
            c.tmp <- -c.tmp
        c.new <- mysd(c.tmp)


        # truncated matrix power update
        ab.tmp <- drop(ttl_rcpp(tnsr = tt, list_mat = list(matrix(c.new, nrow = 1)), ms = c(3)))
        a.new <- mytpower(ab.tmp, sparsity, v1.init = as.vector(a.old), niter = niter)
    }
    return(list(a.new = a.new, c.new = c.new, iter.t = iter.t))
}

#' Function for initial vectors based on aggregation
#' @param T the three order tensor
#' @param sparsity sparsity value at mode 1
#' @param niter maximum number of iterations
#' @keywords internal
myinitial_agg <- function(T, sparsity = NULL) {
    d1 <- dim(T)[1]
    d3 <- dim(T)[3]
    if (is.null(sparsity)) {
        sparsity <- d1
    }
    c0 <- rep(1/sqrt(d3), d3)  # unit initialization
    ab.tmp <- drop(ttl_rcpp(tnsr = T, list_mat = list(matrix(c0, nrow = 1)), ms = c(3)))
    ei <- eigen(ab.tmp)
    ind <- ifelse(abs(ei$values[1]) > abs(ei$values[d1]), 1, d1)
    a.tmp <- ei$vectors[, ind]
    a0 <- mysd(mytruncate(a.tmp, sparsity))
    return(list(a0 = a0))
}

#' Function for initial vectors based on unfolding
#' @param T the three order tensor
#' @param sparsity sparsity value at mode 1
#' @param niter maximum number of iterations
#' @keywords internal
myinitial_unfold <- function(T, sparsity = NULL) {
  d1 <- dim(T)[1]
  d3 <- dim(T)[3]
  if (is.null(sparsity)) {
    sparsity <- d1
  }
  # unfolding
  matT3 <- unfold_tnsr(T,3)
  # the first svd
  aa <- svd(matT3)
  v1 <- aa$v[,1]
  v1m <- matrix(v1,d1,d1)
  # the second svd
  cc <- svd(v1m)
  a0 <- mysd(mytruncate(cc$v[,1], sparsity))
  return(list(a0 = a0))
}



#' Main function for sparse tensor decomposi tion
#' @param T the three order tensor
#' @param sparsity sparsity value at mode 1
#' @param niter maximum number of iterations
#' @param initial the initialization method, can be either aggregation or unfolding.
#' @keywords internal
mystensor <- function(T, sparsity, niter = 10, initial= "aggregation") {
    d1 <- dim(T)[1]
    d3 <- dim(T)[3]

    # initialization
    if(initial == "aggregation"){
      out.initial <- myinitial_agg(T, sparsity = sparsity)
    }else if(initial == "unfolding"){
      out.initial <- myinitial_unfold(T, sparsity = sparsity)
    }


    out.set.all <- myalternative(T, a0 = out.initial$a0, sparsity = sparsity, niter = niter)
    OUT <- list()
    OUT$a <- out.set.all$a.new
    OUT$c <- out.set.all$c.new
    OUT$iter.t <- out.set.all$iter.t
    OUT$initial <- out.initial
    return(OUT)
}



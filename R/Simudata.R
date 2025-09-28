#' Simulate a piecewise stationary  process
#'
#' Simulate piecewise stationary VAR or VMA processes with partial series change and multiple change points. They are used in the simulation studies of Zhang and Chan (2025).
#'
#' @param TT data length.
#' @param m data dimension.
#' @param sparsity number of series activated by change. Cp series are evenly distributed across the `m` dimensions.
#' @param q number of change points, evenly distributed across the `TT` length.
#' @param type VMA or VAR.
#' @return the data `y`.
#' @examples
#' set.seed(1)
#' y <- DGP_mcp(TT = 6000, m = 80, sparsity = 8, q = 4, type = 'vma')
#' @references Zhang and Chan. (2025). Spectral change point estimation for high dimensional time series by sparse tensor decomposition.
#' @export
DGP_mcp <- function(TT, m, sparsity, q, type) {
    cploc <- floor(seq(0, TT, length.out = q + 2))
    seglen <- floor(TT/(q + 1))
    cpseries <- floor(seq(1, m, length.out = sparsity))

    if (type == "vma") {
        sigma <- matrix(0.2, m, m)
        diag(sigma) <- 1
        c1 <- -0.6
        c2 <- -c1
        phi1 <- diag(c1, m)

        phi2 <- phi1
        phi2[cpseries, cpseries] <- diag(c2, sparsity)

        y <- MTS::VARMAsim(nobs = seglen, malags = 1,
            theta = phi1, sigma = sigma)$series
        for (i in 1:q) {
            # the 2nd to (q+1)th segment
            if (i%%2 == 0)
                phi <- phi1 else phi <- phi2
            y2 <- MTS::VARMAsim(nobs = seglen, malags = 1,
                theta = phi, sigma = sigma)$series
            y <- rbind(y, y2)
        }
    } else if (type == "var") {
        sigma <- matrix(0.2, m, m)
        diag(sigma) <- 1
        phi.1 <- diag(0.1, m)

        phi1 <- diag(0.4, m)

        phi2 <- phi1
        phi2[cpseries, cpseries] <- diag(-0.7, sparsity)

        y <- MTS::VARMAsim(nobs = seglen, arlags = c(1, 2), phi = cbind(phi.1, phi1), sigma = sigma)$series  #note arlags
        for (i in 1:q) {
            # the 2nd to (q+1)th segment
            if (i%%2 == 0)
                phi <- phi1 else phi <- phi2
            y2 <- MTS::VARMAsim(nobs = seglen, arlags = c(1,
                2), phi = cbind(phi.1, phi), sigma = sigma)$series
            y <- rbind(y, y2)
        }
    }
    return(y)
}

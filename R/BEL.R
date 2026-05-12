#' The two-sample blockwise empirical likelihood statistic for differences in means
#'
#' @description Calculates blockwise empirical likelihood test for the difference of two sample means.
#'
#' @param X,Y vectors of data values.
#' @param M_1,M_2 positive integers specifying block length for X and Y, respectively.
#' @param Delta hypothesized difference of two populations.
#'
#' @return A list of class "htest" containing following components:
#' method - the character string of the test.
#' data.name - a character string with the names of the input data.
#' Delta0 - the specified hypothesized value of mean differences under the null hypothesis
#' statistic - the value of the test statistic.
#' p.value - the p-value for the test.
#'
#' @author R. Alksnis, J. Valeinis
#'
#' @examples
#' # Basic example
#' Delta0 <- 1.5
#' X <- arima.sim(n = 400, model = list(ar = .3))
#' Y <- arima.sim(n = 400, model = list(ar = .5)) + Delta0
#' BEL.means(X, Y, 10, 20, Delta = Delta0)
#'
#' @export
BEL.means <- function(X, Y, M_1, M_2, Delta = 0) {
    blocks <- function(x, M) {
        if (!((M %% 1 == 0) * (M > 0))) {
            stop("Block lengths M_1 and M_2 need to be positive integers!")
        }
        N <- length(x)
        Q <- floor((N - M) / M) + 1
        means <- c()
        for (i in 1:Q) {
            means[i] <- mean(x[((i - 1) * M + 1):(M * i)])
        }
        return(means)
    }

    system <- function(x, Delta, sampleX, sampleY, M_1, M_2) {
        y <- numeric(3)
        N_1 <- length(sampleX)
        N_2 <- length(sampleY)
        T_i <- blocks(sampleX - x[3], M_1)
        T_j <- blocks(sampleY - x[3] - Delta, M_2)
        y[1] <- sum(T_i / (1 + x[1] * T_i))
        y[2] <- sum(T_j / (1 + x[2] * T_j))
        y[3] <- x[1] * sum(1 / (1 + x[1] * T_i)) + x[2] * sum(1 / (1 + x[2] * T_j))
        return(y)
    }
    N_1 <- length(X)
    N_2 <- length(Y)
    Q_1 <- floor((N_1 - M_1) / M_1) + 1
    Q_2 <- floor((N_2 - M_2) / M_2) + 1
    x_start <- c(0, 0, mean(X))
    #' @import nleqslv
    solution <- nleqslv::nleqslv(
        x_start,
        system,
        Delta = Delta,
        sampleX = X,
        sampleY = Y,
        M_1 = M_1,
        M_2 = M_2
    )
    lambda_1_star <- solution$x[1]
    lambda_2_star <- solution$x[2]
    mu_star <- solution$x[3]
    Tx_i <- blocks(X - mu_star, M_1)
    Ty_j <- blocks(Y - mu_star - Delta, M_2)
    ELR <-  (2 * sum(log(1 + lambda_1_star * Tx_i)) + 2 * sum(log(1 + lambda_2_star * Ty_j)))
    p <- 1 - pchisq(ELR, df = 1)
    test <- "Two-Sample Blockwise Empirical Likelihood (BEL) test"
    attr(ELR, 'names') <- '-2 * LogLikelihood'
    attr(p, 'names') <- NULL
    attr(Delta, "names") <- "Difference of means"
    data  <- paste0(deparse(substitute(X)), " and ", deparse(substitute(Y)))
    alternative <- 'two.sided'
    Results  <- list(
        method = test,
        data.name = data,
        null.value = Delta,
        alternative = alternative,
        statistic = ELR,
        p.value = p
    )
    class(Results) <- "htest"
    Results
}

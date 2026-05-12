#' Two-sample frequency-domain empirical likelihood test for
#' autocorrelation difference
#'
#' Tests whether the lag-\eqn{h} autocorrelation of two
#' independent stationary time series differ by a specified
#' amount \eqn{\Delta}.
#'
#' @param X,Y Numeric time-series vectors (length >= 10).
#' @param Delta Hypothesised difference
#'   \eqn{\rho_x(h) - \rho_y(h) = \Delta} (default 0).
#' @param lag Positive integer lag \eqn{h} (default 1).
#' @param bartlett Logical; apply a bootstrap Bartlett
#'   correction? (default \code{FALSE}).
#' @param bootstrap.samples Number of bootstrap replicates for the Bartlett
#'   correction (used only when \code{bartlett = TRUE}).
#' @param center Logical; subtract sample means before
#'   computing periodograms? (default \code{TRUE}).
#' @param seed Optional integer seed for reproducibility.
#' @param span Loess span for periodogram smoothing used in the
#'   Bartlett correction (default 0.15).
#' @param rho.lower,rho.upper Search bounds for the profile
#'   autocorrelation optimization (defaults \code{-0.99},
#'   \code{0.99}).
#'
#' @return An object of class \code{c("FDELacf", "htest")}
#'   with components \code{statistic}, \code{parameter},
#'   \code{p.value}, \code{estimate}, \code{null.value},
#'   \code{alternative}, \code{method}, \code{data.name},
#'   \code{lag}, \code{bartlett}, \code{bartlett.factor},
#'   \code{statistic.uncorrected}, \code{p.value.uncorrected},
#'   \code{B}, \code{call}.
#'
#' @author R. Alksnis, J. Valeinis
#'
#' @examples
#' set.seed(1)
#' X <- arima.sim(n = 200, model = list(ar = 0.3))
#' Y <- arima.sim(n = 350, model = list(ar = 0.3))
#' FDEL.acf(X, Y)
#' ## With Bartlett correction (slower)
#' ## FDEL.acf(X, Y, bartlett = TRUE, B = 199)
#'
#' @importFrom stats pchisq loess predict uniroot optimize fft median rexp runif
#' @export
FDEL.acf <- function(X,
                     Y,
                     Delta = 0,
                     lag = 1,
                     bartlett = FALSE,
                     bootstrap.samples = 500,
                     center = TRUE,
                     seed = NULL,
                     span = 0.15,
                     rho.lower = -0.99,
                     rho.upper = 0.99) {
    B <- bootstrap.samples

    if (!is.numeric(X) || !is.numeric(Y)) {
        stop("'X' and 'Y' must be numeric.")
    }

    if (length(X) < 10 || length(Y) < 10) {
        stop("Series are too short.")
    }

    if (!is.numeric(Delta) ||
        length(Delta) != 1 || !is.finite(Delta)) {
        stop("'Delta' must be a finite numeric value.")
    }

    if (!is.numeric(lag) ||
        length(lag) != 1 || lag < 1 || lag != floor(lag)) {
        stop("'lag' must be a positive integer.")
    }

    if (!is.logical(bartlett) || length(bartlett) != 1) {
        stop("'bartlett' must be TRUE or FALSE.")
    }

    if (bartlett &&
        (!is.numeric(B) ||
         length(B) != 1 || B < 1 || B != floor(B))) {
        stop("'bootstrap.samples' must be a positive integer when bartlett = TRUE.")
    }

    if (!is.numeric(span) ||
        length(span) != 1 || span <= 0 || span > 1) {
        stop("'span' must be in (0, 1].")
    }

    if (!is.numeric(rho.lower) || !is.numeric(rho.upper) ||
        length(rho.lower) != 1 || length(rho.upper) != 1 ||
        !is.finite(rho.lower) || !is.finite(rho.upper) ||
        rho.lower >= rho.upper) {
        stop("'rho.lower' and 'rho.upper' must be finite with rho.lower < rho.upper.")
    }

    if (!is.null(seed)) {
        had_seed <- exists(".Random.seed",
                           envir = .GlobalEnv,
                           inherits = FALSE)

        if (had_seed) {
            old_seed <- get(".Random.seed", envir = .GlobalEnv)
        }

        on.exit({
            if (had_seed) {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
            } else if (exists(".Random.seed",
                              envir = .GlobalEnv,
                              inherits = FALSE)) {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        }, add = TRUE)

        set.seed(seed)
    }

    X <- as.numeric(X)
    Y <- as.numeric(Y)

    if (anyNA(X) || anyNA(Y)) {
        stop("'X' and 'Y' must not contain missing values.")
    }

    if (!all(is.finite(X)) || !all(is.finite(Y))) {
        stop("'X' and 'Y' must contain only finite values.")
    }

    if (center) {
        X <- X - mean(X)
        Y <- Y - mean(Y)
    }

    stat_raw <- .fdel_acf_stat(
        X,
        Y,
        Delta = Delta,
        lag = lag,
        center = FALSE,
        rho.lower = rho.lower,
        rho.upper = rho.upper
    )

    ell_uncorrected <- stat_raw$statistic

    if (!is.finite(ell_uncorrected)) {
        stop(
            "FDEL statistic could not be computed. The empirical likelihood constraints may be infeasible."
        )
    }

    pval_uncorrected <- 1 - stats::pchisq(ell_uncorrected, df = 1)

    ell <- ell_uncorrected
    pval <- pval_uncorrected
    tau <- NA_real_

    if (bartlett) {
        tau <- .fdel_acf_bartlett_factor(
            X,
            Y,
            Delta = Delta,
            lag = lag,
            rho0 = stat_raw$rho.hat,
            B = B,
            span = span,
            rho.lower = rho.lower,
            rho.upper = rho.upper
        )
        ell  <- ell_uncorrected / tau
        pval <- 1 - stats::pchisq(ell, df = 1)
    }

    out <- list(
        statistic   = c("FDEL chi-square" = ell),
        parameter   = c(df = 1),
        p.value     = pval,
        estimate    = c("profile rho" = stat_raw$rho.hat),
        null.value  = c(Delta = Delta),
        alternative = "two.sided",
        method      = if (bartlett)
            paste(
                "Two-sample FDEL test for autocorrelation difference",
                "with bootstrap Bartlett correction"
            )
        else
            "Two-sample FDEL test for autocorrelation difference",
        data.name              = paste(deparse(substitute(X)),
                                       "and",
                                       deparse(substitute(Y))),
        lag                    = lag,
        bartlett               = bartlett,
        bartlett.factor        = tau,
        statistic.uncorrected  = c("FDEL chi-square" = ell_uncorrected),
        p.value.uncorrected    = pval_uncorrected,
        bootstrap.samples      = if (bartlett) B else NA_integer_,
        call                   = match.call()
    )
    class(out) <- c("FDELacf", "htest")
    out
}

#' @export
print.FDELacf <- function(x, digits = getOption("digits"), ...) {
    y <- x
    class(y) <- "htest"
    print(y, digits = digits, ...)
    cat("\nLag:", x$lag, "\n")
    if (isTRUE(x$bartlett)) {
        cat("\nBootstrap Bartlett correction:\n")
        cat("Bartlett factor:        ",
            format(x$bartlett.factor, digits = digits),
            "\n",
            sep = "")
        cat("Bootstrap replications: ",
            x$bootstrap.samples,
            "\n",
            sep = "")
        cat("Uncorrected statistic:  ",
            format(unname(x$statistic.uncorrected), digits = digits),
            "\n",
            sep = "")
        cat(
            "Uncorrected p-value:    ",
            format.pval(x$p.value.uncorrected, digits = digits),
            "\n",
            sep = ""
        )
    }
    invisible(x)
}

.fdel_acf_stat <- function(X,
                           Y,
                           Delta = 0,
                           lag = 1,
                           center = TRUE,
                           rho.lower = -0.99,
                           rho.upper = 0.99) {
    if (center) {
        X <- X - mean(X)
        Y <- Y - mean(Y)
    }
    n <- length(X)
    m <- length(Y)
    N <- floor((n - 1) / 2)
    M <- floor((m - 1) / 2)

    IX <- abs(fft(X))^2 / (2 * pi * n)
    IY <- abs(fft(Y))^2 / (2 * pi * m)

    lX <- 2 * pi * seq_len(N) / n
    lY <- 2 * pi * seq_len(M) / m

    IX <- IX[2:(N + 1)]
    IY <- IY[2:(M + 1)]

    objective <- function(rho) {
        gX <- (cos(lag * lX) - rho) * IX
        gY <- (cos(lag * lY) - rho - Delta) * IY
        lamX <- .fdel_solve_lambda(gX)
        lamY <- .fdel_solve_lambda(gY)
        if (is.null(lamX) || is.null(lamY))
            return(.Machine$double.xmax)
        val <- 4 * (sum(log1p(pi * lamX * gX)) +
                        sum(log1p(pi * lamY * gY)))
        if (!is.finite(val))
            .Machine$double.xmax
        else
            val
    }

    opt <- optimize(f = objective,
                    lower = rho.lower,
                    upper = rho.upper)
    list(statistic = opt$objective,
         rho.hat = opt$minimum)
}

.fdel_solve_lambda <- function(g) {
    if (!all(is.finite(g)))
        return(NULL)
    if (!(min(g) < 0 && max(g) > 0))
        return(NULL)
    u <- pi * g
    lower <- if (any(u > 0))
        max(-1 / u[u > 0]) + 1e-8
    else
        - 1e8
    upper <- if (any(u < 0))
        min(-1 / u[u < 0]) - 1e-8
    else
        1e8
    if (!(lower < upper))
        return(NULL)
    f <- function(lambda) {
        denom <- 1 + lambda * u
        if (any(denom <= 0))
            return(NA_real_)
        sum(g / denom)
    }
    sol <- tryCatch(
        uniroot(
            f,
            lower = lower,
            upper = upper,
            tol = .Machine$double.eps^0.25
        ),
        error = function(e)
            NULL
    )
    if (is.null(sol))
        return(NULL)
    lambda <- sol$root
    if (any(1 + lambda * u <= 0))
        return(NULL)
    lambda
}

.fdel_acf_bartlett_factor <- function(X,
                                      Y,
                                      Delta = 0,
                                      lag = 1,
                                      rho0,
                                      B = 500,
                                      span = 0.15,
                                      rho.lower = -0.99,
                                      rho.upper = 0.99) {
    n <- length(X)
    m <- length(Y)
    N <- floor((n - 1) / 2)
    M <- floor((m - 1) / 2)

    IX <- abs(fft(X))^2 / (2 * pi * n)
    IY <- abs(fft(Y))^2 / (2 * pi * m)
    IX <- IX[2:(N + 1)]
    IY <- IY[2:(M + 1)]

    lX <- 2 * pi * seq_len(N) / n
    lY <- 2 * pi * seq_len(M) / m

    fX <- .fdel_smooth_periodogram(IX, span = span)
    fY <- .fdel_smooth_periodogram(IY, span = span)

    fX0 <- .fdel_tilt_spectrum(lX, fX, rho0, lag)
    fY0 <- .fdel_tilt_spectrum(lY, fY, rho0 + Delta, lag)

    ell.star <- numeric(B)
    for (b in seq_len(B)) {
        Xs <- .fdel_pseudo_series(fX0, n)
        Ys <- .fdel_pseudo_series(fY0, m)
        ell.star[b] <- tryCatch(
            .fdel_acf_stat(
                Xs,
                Ys,
                Delta = Delta,
                lag = lag,
                center = TRUE,
                rho.lower = rho.lower,
                rho.upper = rho.upper
            )$statistic,
            error = function(e)
                NA_real_
        )
    }
    ell.star[ell.star >= .Machine$double.xmax / 10] <- NA_real_
    tau <- mean(ell.star, na.rm = TRUE)
    if (!is.finite(tau) || tau <= 0) {
        warning("Bootstrap Bartlett factor could not be estimated; using factor 1.")
        tau <- 1
    }
    tau
}

.fdel_smooth_periodogram <- function(I, span = 0.15) {
    I <- pmax(I, 1e-12)
    k <- seq_along(I)
    fit <- tryCatch(
        stats::loess(
            log(I) ~ k,
            span = span,
            degree = 1,
            surface = "direct"
        ),
        error = function(e)
            NULL
    )
    if (is.null(fit))
        return(rep(mean(I), length(I)))
    out <- exp(predict(fit, newdata = data.frame(k = k)))
    out[!is.finite(out)] <- median(I)
    pmax(out, 1e-12)
}

.fdel_tilt_spectrum <- function(freq, fhat, target.rho, lag = 1) {
    target.rho <- max(min(target.rho, 0.999), -0.999)
    moment <- function(theta) {
        w <- fhat * exp(theta * cos(lag * freq))
        sum(cos(lag * freq) * w) / sum(w) - target.rho
    }
    lower <- -100
    upper <- 100
    ml <- moment(lower)
    mu <- moment(upper)
    if (!is.finite(ml) || !is.finite(mu) || ml * mu > 0)
        return(fhat)
    theta <- uniroot(moment, lower = lower, upper = upper)$root
    pmax(fhat * exp(theta * cos(lag * freq)), 1e-12)
}

.fdel_pseudo_series <- function(f0, n) {
    N    <- floor((n - 1) / 2)
    E    <- rexp(N)
    Istar <- f0 * E
    amp   <- sqrt(2 * pi * n * Istar)
    phase <- runif(N, 0, 2 * pi)
    coef.pos <- amp * exp(1i * phase)
    Z <- complex(length.out = n)
    Z[1] <- 0
    Z[2:(N + 1)] <- coef.pos
    Z[(n - N + 1):n] <- Conj(rev(coef.pos))
    if (n %% 2 == 0)
        Z[n / 2 + 1] <- Re(Z[n / 2 + 1])
    Re(fft(Z, inverse = TRUE) / n)
}

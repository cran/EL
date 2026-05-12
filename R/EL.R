#' The two-sample empirical likelihood statistic
#'
#' Calculates \eqn{-2} times the log-likelihood ratio statistic when the function
#' of interest (either of P-P or Q-Q plot, ROC curve, difference of quantile or
#' distribution functions) at some point \code{t} is equal to \code{d}.
#'
#' @param method One of \code{"pp"}, \code{"qq"}, \code{"roc"}, \code{"qdiff"},
#'   or \code{"fdiff"}.
#' @param X A numeric vector of data values.
#' @param Y A numeric vector of data values.
#' @param Delta A number. The hypothesised value of the function at point \code{t}.
#' @param d Deprecated. Use \code{Delta} instead.
#' @param t A number. The point at which to evaluate the statistic.
#' @param bw A function taking a numeric vector and returning a bandwidth, or a
#'   numeric vector of length two giving the bandwidths for \code{X} and \code{Y}
#'   respectively. Default is \code{\link[stats]{bw.nrd0}}.
#' @param conf.level Confidence level for the reported
#'   confidence interval (default \code{0.95}).
#'
#' @return An object of class \code{"htest"}.
#'
#' @references
#' Valeinis, J. and Cers, E. Extending the two-sample empirical likelihood.
#' Preprint: \url{http://home.lu.lv/~valeinis/lv/petnieciba/EL_TwoSample_2011.pdf}
#'
#' @author E. Cers, J. Valeinis
#'
#' @seealso \code{\link{EL.smooth}}
#'
#' @examples
#' EL.statistic(method = "pp", X = rnorm(100), Y = rnorm(100), Delta = 0.5, t = 0.5)
#'
#' @keywords nonparametric smooth htest
#' @importFrom stats bw.nrd0 dnorm ecdf optimize pchisq pnorm qchisq qnorm quantile uniroot
#' @export
EL.statistic <- function(method,
                         X,
                         Y,
                         Delta,
                         d,
                         t,
                         bw = bw.nrd0,
                         conf.level = 0.95) {
    if (!missing(d)) {
        warning("Argument 'd' is deprecated, use 'Delta' instead.",
                call. = FALSE)
        if (missing(Delta))
            Delta <- d
    }
    cl <- match.call()

    stat <- EL.local(ELLR, X, Y, Delta, t, method = method, bw = bw)

    ci <- tryCatch(
        EL.local(
            function(globals, X, Y)
                ELconf(globals, X, Y, t, conf.level),
            X,
            Y,
            method = method,
            bw = bw
        ),
        error = function(e)
            c(NA_real_, NA_real_)
    )
    attr(ci, "conf.level") <- conf.level

    method_str <- switch(
        method,
        pp    = "EL test for P-P plot",
        qq    = "EL test for Q-Q plot",
        roc   = "EL test for ROC curve",
        qdiff = "EL test for quantile difference",
        fdiff = "EL test for CDF difference",
        paste("EL test -", method)
    )

    structure(
        list(
            statistic  = c("-2 log LR" = stat),
            parameter  = c(df = 1),
            p.value    = 1 - pchisq(stat, df = 1),
            estimate   = c(estimate = Delta),
            conf.int   = ci,
            null.value = c(parameter = Delta),
            alternative = "two.sided",
            method     = method_str,
            data.name  = paste(deparse(cl$X), "and", deparse(cl$Y))
        ),
        class = "htest"
    )
}

#' Draws plots using the smoothed two-sample empirical likelihood method
#'
#' Draws P-P and Q-Q plots, ROC curves, quantile differences (\code{qdiff}) and
#' CDF differences (\code{ddiff}) and their respective confidence bands (pointwise
#' or simultaneous) using the empirical likelihood method.
#'
#' @inheritParams EL.smooth
#' @param ... Further arguments passed to the plot function.
#'
#' @details
#' The plotting interval for P-P plots, ROC curves and differences of quantile
#' functions is \eqn{[0, 1]}. The Q-Q plot is drawn from the minimum to the
#' maximum of \code{Y}. For the plot of distribution function differences the
#' interval from \eqn{\max(\min(X), \min(Y))} to \eqn{\min(\max(X), \max(Y))}
#' is used.
#'
#' Confidence bands are drawn only if \code{conf.level} is not \code{NULL}.
#'
#' When constructing simultaneous confidence bands, the plot is drawn on an
#' interval narrowed by 5\% on both sides, since the procedure is usually
#' sensitive at the endpoints. The confidence level is bootstrapped using 50
#' evenly spaced points in this interval. If the default interval produces
#' excessively wide bands, use \code{\link{EL.smooth}} where intervals are
#' specified manually. Note that calculation of simultaneous confidence bands
#' can be slow.
#'
#' @return A \code{ggplot} object. The plot can be further customized using
#'   \pkg{ggplot2} functions, as shown in the examples.
#'
#' @references
#' Valeinis, J. and Cers, E. Extending the two-sample empirical likelihood.
#' Preprint: \url{http://home.lu.lv/~valeinis/lv/petnieciba/EL_TwoSample_2011.pdf}
#'
#' Hall, P. and Owen, A. (1993). Empirical likelihood bands in density estimation.
#' \emph{Journal of Computational and Graphical Statistics}, \strong{2}(3), 273--289.
#'
#' @author E. Cers, J. Valeinis
#'
#' @seealso \code{\link{EL.smooth}}, \code{\link{EL.statistic}}
#'
#' @examples
#' ## The examples showcase all available graphs
#' set.seed(42)
#' X1 <- rnorm(100, 0.5, 1.5)
#' X2 <- rnorm(100, 0, 1)
#'
#' \donttest{
#' xlim <- c(min(X1, X2) - 0.5, max(X1, X2) + 0.5)
#' D1 <- density(X1)
#' D2 <- density(X2)
#' df <- data.frame(x1 = D1$x, y1 = D1$y, x2 = D2$x, y2 = D2$y)
#' p1 <- ggplot2::ggplot(data = df) +
#'     ggplot2::geom_line(ggplot2::aes(x = x2, y = y2,
#'         color = paste0('X2 (bw=', round(D2$bw, 2), ')'))) +
#'     ggplot2::geom_line(ggplot2::aes(x = x1, y = y1,
#'         color = paste0('X1 (bw=', round(D1$bw, 2), ')'))) +
#'     ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) +
#'     ggplot2::theme_minimal() +
#'     ggplot2::theme(legend.position = "top") +
#'     ggplot2::labs(x = "X", y = "Density")
#' p1
#'
#' # CDF differences
#' p2 <- EL.plot("fdiff", X1, X2, main = "F difference", conf.level = 0.95)
#' tt <- seq(max(c(min(X1), min(X2))), min(c(max(X1), max(X2))), length = 30)
#' ee <- ecdf(X2)(tt) - ecdf(X1)(tt)
#' p2 <- p2 + ggplot2::geom_point(data = data.frame(tt = tt, ee = ee),
#'                                 ggplot2::aes(x = tt, y = ee))
#' p2
#'
#' # Quantile differences
#' p3 <- EL.plot("qdiff", X1, X2, main = "Quantile difference", conf.level = 0.95)
#' tt <- seq(0.01, 0.99, length = 30)
#' ee <- quantile(X2, tt) - quantile(X1, tt)
#' p3 <- p3 + ggplot2::geom_point(data = data.frame(tt = tt, ee = ee),
#'                                 ggplot2::aes(x = tt, y = ee))
#' p3
#'
#' # Q-Q plot
#' p4 <- EL.plot("qq", X1, X2, main = "Q-Q plot", conf.level = 0.95)
#' tt <- seq(min(X2), max(X2), length = 30)
#' ee <- quantile(X1, ecdf(X2)(tt))
#' p4 <- p4 + ggplot2::geom_point(data = data.frame(tt = tt, ee = ee),
#'                                 ggplot2::aes(x = tt, y = ee))
#' p4
#'
#' # P-P plot
#' p5 <- EL.plot("pp", X1, X2, main = "P-P plot", conf.level = 0.95, ylim = c(0, 1))
#' tt <- seq(0.01, 0.99, length = 30)
#' ee <- ecdf(X1)(quantile(X2, tt))
#' p5 <- p5 + ggplot2::geom_point(data = data.frame(tt = tt, ee = ee),
#'                                 ggplot2::aes(x = tt, y = ee))
#' p5
#'
#' # ROC curve
#' p6 <- EL.plot("roc", X1, X2, main = "ROC curve", conf.level = 0.95, ylim = c(0, 1))
#' tt <- seq(0.01, 0.99, length = 30)
#' ee <- 1 - ecdf(X1)(quantile(X2, 1 - tt))
#' p6 <- p6 + ggplot2::geom_point(data = data.frame(tt = tt, ee = ee),
#'                                 ggplot2::aes(x = tt, y = ee))
#' p6
#'
#' # To show all plots at once:
#' # require(cowplot)
#' # cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
#' }
#'
#' @keywords nonparametric smooth hplot
#' @export
EL.plot <- function(method,
                    X,
                    Y,
                    bw = bw.nrd0,
                    conf.level = NULL,
                    simultaneous = FALSE,
                    bootstrap.samples = 300,
                    more.warnings = FALSE,
                    ...)
{
    numpoints <- 50
    call <- match.call(expand.dots = FALSE)
    optional <- call$...
    switch(match(method, c("fdiff", "qdiff", "qq", "pp", "roc")), {
        dxlab <- "x"
        dylab <- substitute(expression(F[a] - F[b]), list(a = call$Y, b = call$X))
    }, {
        dxlab <- "P"
        dylab <- substitute(expression(F[a]^-1 - F[b]^-1), list(a = call$Y, b = call$X))
    }, {
        dxlab <- substitute(expression(F[a]^-1), list(a = call$Y))
        dylab <- substitute(expression(F[a]^-1), list(a = call$X))
    }, {
        dxlab <- substitute(expression(F[a]), list(a = call$Y))
        dylab <- substitute(expression(F[a]), list(a = call$X))
    }, {
        dxlab <- "False positive rate"
        dylab <- "True positive rate"
    })
    tr <- EL.local(function(globals, X, Y)
        globals$trange(globals, X, Y),
        X,
        Y,
        method = method,
        bw = bw)
    if (simultaneous)
    {
        trlen <- tr[2] - tr[1]
        tr[1] <- tr[1] + 0.05 * trlen
        tr[2] <- tr[2] - 0.05 * trlen
    }
    t <- seq(tr[1], tr[2], length = numpoints)
    zz <- EL.smooth(
        method = method,
        X,
        Y,
        t,
        bw = bw,
        conf.level = conf.level,
        simultaneous = simultaneous,
        bootstrap.samples = bootstrap.samples,
        more.warnings = more.warnings
    )
    if (!is.null(conf.level))
        dylim <- c(min(zz$conf.int[1, ], na.rm = TRUE),
                   max(zz$conf.int[2, ], na.rm = TRUE))
    else
        dylim <- c(min(zz$estim), max(zz$estim))

    # Avoid checking issues due to ggplot calls
    estim <- lower_conf <- upper_conf <- NULL

    # Create a data frame for ggplot
    df <- data.frame(
        t = t,
        estim = zz$estim,
        lower_conf = if (!is.null(conf.level))
            zz$conf.int[1, ]
        else
            NA,
        upper_conf = if (!is.null(conf.level))
            zz$conf.int[2, ]
        else
            NA
    )

    xlab_val <- if ("xlab" %in% names(optional)) {
        optional$xlab
    } else if (is.character(dxlab)) {
        dxlab
    } else {
        eval(dxlab)
    }

    ylab_val <- if ("ylab" %in% names(optional)) {
        optional$ylab
    } else if (is.character(dylab)) {
        dylab
    } else {
        eval(dylab)
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = t, y = estim)) +
        ggplot2::geom_line(color = "blue") +
        ggplot2::labs(x = xlab_val, y = ylab_val) +
        ggplot2::theme_minimal()

    if (!is.null(conf.level)) {
        p <- p +
            ggplot2::geom_ribbon(
                ggplot2::aes(ymin = lower_conf, ymax = upper_conf),
                fill = "grey80",
                alpha = 0.5
            ) +
            ggplot2::geom_line(ggplot2::aes(y = lower_conf),
                               linetype = "dashed",
                               color = "red") +
            ggplot2::geom_line(ggplot2::aes(y = upper_conf),
                               linetype = "dashed",
                               color = "red")
    }
    p
}

#' Smooth estimates and confidence intervals using the smoothed two-sample EL method
#'
#' Calculates estimates and pointwise confidence intervals (or simultaneous bands)
#' for P-P and Q-Q plots, ROC curves, quantile differences (\code{qdiff}) and CDF
#' differences (\code{ddiff}) using the smoothed empirical likelihood method.
#'
#' @inheritParams EL.statistic
#' @param t A numeric vector of points at which to calculate estimates and
#'   confidence intervals.
#' @param conf.level Confidence level for the intervals, a number in \eqn{(0,1)},
#'   or \code{NULL} to skip confidence band calculation.
#' @param simultaneous Logical. If \code{TRUE}, simultaneous confidence bands are
#'   constructed via a nonparametric bootstrap. Default is \code{FALSE} (pointwise
#'   intervals).
#' @param bootstrap.samples Integer. Number of bootstrap samples used when
#'   \code{simultaneous = TRUE}. Default is \code{300}.
#' @param more.warnings Logical. If \code{FALSE} (default), a single warning is
#'   issued if any problem occurs. If \code{TRUE}, a warning is produced for every
#'   point with a problem.
#'
#' @details
#' Confidence bands are drawn only if \code{conf.level} is not \code{NULL}.
#'
#' When constructing simultaneous confidence bands, check that the chosen range
#' of \code{t} values does not produce excessively wide bands (for example, for
#' a P-P plot the interval \eqn{[0.05, 0.95]} is typically a sensible choice).
#' This should be verified for each dataset. Note that simultaneous band
#' calculation can be slow.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{estimate}}{Estimated values at points \code{t}.}
#'   \item{\code{conf.int}}{A two-row matrix where each column gives the lower
#'     and upper confidence bounds at the corresponding point in \code{t}.}
#'   \item{\code{simultaneous.conf.int}}{Logical; \code{TRUE} if simultaneous
#'     bands were constructed.}
#'   \item{\code{bootstrap.crit}}{The bootstrap critical value of the
#'     \eqn{-2\log}-likelihood statistic for simultaneous bands at level
#'     \code{conf.level}. Only present when \code{conf.level} is not \code{NULL}
#'     and \code{simultaneous = TRUE}.}
#' }
#'
#' @references
#' Valeinis, J. and Cers, E. Extending the two-sample empirical likelihood.
#' Preprint: \url{http://home.lu.lv/~valeinis/lv/petnieciba/EL_TwoSample_2011.pdf}
#'
#' Hall, P. and Owen, A. (1993). Empirical likelihood bands in density estimation.
#' \emph{Journal of Computational and Graphical Statistics}, \strong{2}(3), 273--289.
#'
#' @author E. Cers, J. Valeinis
#'
#' @seealso \code{\link{EL.plot}}, \code{\link{EL.statistic}}
#'
#' @examples
#' \donttest{
#' #### Simultaneous confidence bands for a P-P plot
#' X1 <- rnorm(200)
#' X2 <- rnorm(200, 1)
#' x <- seq(0.05, 0.95, length = 19)
#' y <- EL.smooth("pp", X1, X2, x, conf.level = 0.95,
#'                simultaneous = TRUE, bw = c(0.3, 0.3))
#' conf.int <- data.frame(x = x, ci.l = y$conf.int[1,], ci.u = y$conf.int[2,])
#'
#' ## Plot with both pointwise and simultaneous confidence bands
#' EL.plot("pp", X1, X2, conf.level = 0.95, bw = c(0.3, 0.3)) +
#'     ggplot2::geom_line(data = conf.int, ggplot2::aes(x = x, y = ci.u), lty = "dotted") +
#'     ggplot2::geom_line(data = conf.int, ggplot2::aes(x = x, y = ci.l), lty = "dotted")
#' }
#'
#' @keywords hplot nonparametric smooth
#' @export
EL.smooth <- function(method,
                      X,
                      Y,
                      t,
                      bw = bw.nrd0,
                      conf.level = NULL,
                      simultaneous = FALSE,
                      bootstrap.samples = 300,
                      more.warnings = FALSE)
{
    warnings <- FALSE
    mm <- function(globals, tt)
    {
        d <- sapply(tt, function(t)
            deltasolve(globals, X, Y, t))
        bs <- NA
        conffun <- function(conf.level)
        {
            function(k)
                withCallingHandlers({
                    ci <- ELconf(globals, X, Y, tt[k], conf.level, delta = d[k])
                    if (!any(is.na(ci))) {
                        ci[1] <- min(ci[1], d[k])
                        ci[2] <- max(ci[2], d[k])
                    }
                    ci
                }, simpleWarning = function(e)
                {
                    if (!more.warnings)
                    {
                        warnings <<- TRUE
                        invokeRestart("muffleWarning")
                    }
                })
        }
        if (is.null(conf.level))
            ci <- NA
        else
        {
            if (!simultaneous)
            {
                ci <- sapply(seq_along(tt), conffun(conf.level))
                attr(ci, "conf.level") <- conf.level
            }
            else
            {
                bootup <- function()
                {
                    X1 <- sample(X, replace = TRUE)
                    Y1 <- sample(Y, replace = TRUE)
                    vals <- max(mapply(
                        ELLR,
                        delta = d,
                        t = tt,
                        MoreArgs = list(
                            globals = globals,
                            X = X1,
                            Y = Y1
                        )
                    ))
                    vals[!is.finite(vals)] <- NA_real_
                    max(vals, na.rm = TRUE)
                }
                bs <- sort(replicate(bootstrap.samples, bootup()))
                bs <- bs[ceiling(conf.level * bootstrap.samples)]
                ci <- sapply(seq_along(tt), conffun(pchisq(bs, 1)))
                for (k in seq_along(tt)) {
                    if (!any(is.na(ci[, k]))) {
                        ci[1, k] <- min(ci[1, k], d[k])
                        ci[2, k] <- max(ci[2, k], d[k])
                    }
                }
            }
        }

        list(
            estimate = d,
            conf.int = ci,
            simultaneous.conf.int = simultaneous,
            bootstrap.crit = bs
        )
    }
    res <- EL.local(mm, t, method = method, bw = bw)
    if (warnings)
        warning(
            paste(
                "Estimates or confidence bands for some values could",
                "not be found. For a detailed report, run this command",
                "again with 'more.warnings = TRUE'."
            ),
            call. = FALSE
        )
    res
}

#' Empirical likelihood test for the difference of two sample means
#'
#' Empirical likelihood inference for the difference of two sample means.
#' This includes a test for the null hypothesis of a constant mean difference,
#' a confidence interval, and the EL estimator.
#'
#' @inheritParams EL.statistic
#' @param mu A number specifying the null hypothesis value for the mean
#'   difference. Default is \code{0}.
#'
#' @return A list of class \code{"htest"} with components:
#' \describe{
#'   \item{\code{estimate}}{The empirical likelihood estimate of the mean difference.}
#'   \item{\code{conf.int}}{A confidence interval for the mean difference.}
#'   \item{\code{p.value}}{The p-value for the test.}
#'   \item{\code{statistic}}{The value of the test statistic.}
#'   \item{\code{method}}{The character string \code{"Empirical likelihood mean difference test"}.}
#'   \item{\code{null.value}}{The hypothesised mean difference \code{mu} under \eqn{H_0}.}
#'   \item{\code{data.name}}{A character string giving the names of the data.}
#' }
#'
#' @references
#' Valeinis, J., Cers, E. and Cielens, J. (2010). Two-sample problems in statistical
#' data modelling. \emph{Mathematical Modelling and Analysis}, \strong{15}(1), 137--151.
#'
#' Valeinis, J. and Cers, E. Extending the two-sample empirical likelihood.
#' Preprint: \url{http://home.lu.lv/~valeinis/lv/petnieciba/EL_TwoSample_2011.pdf}
#'
#' @author E. Cers, J. Valeinis
#'
#' @seealso \code{\link{EL.Huber}}
#'
#' @examples
#' X <- rnorm(100)
#' Y <- rnorm(100)
#' t.test(X, Y)
#' EL.means(X, Y)
#' EL.Huber(X, Y)
#'
#' @keywords nonparametric htest
#' @export
EL.means <- function(X,
                     Y,
                     mu = 0,
                     conf.level = 0.95) {
    call <- match.call()
    mm <- function(globals)
    {
        d <- deltasolve(globals, X, Y, 0)
        names(d) <- "Mean difference"
        ci <- ELconf(globals, X, Y, 0, conf.level, delta = d)
        attr(ci, "conf.level") <- conf.level
        ll <- ELLR(globals, X, Y, mu, 0)
        names(ll) <- "-2 log LR"
        res <- list(
            statistic   = ll,
            parameter   = c(df = 1),
            p.value     = 1 - pchisq(ll, 1),
            estimate    = d,
            conf.int    = ci,
            null.value  = c("mean difference" = mu),
            alternative = "two.sided",
            method      = "Two-sample empirical likelihood mean difference test",
            data.name   = paste(deparse(call$X), "and", deparse(call$Y))
        )
        class(res) <- "htest"
        res
    }
    EL.local(mm, method = "mean")
}

#' Empirical likelihood test for the difference of smoothed Huber estimators
#'
#' Empirical likelihood inference for the difference of smoothed Huber estimators.
#' This includes a test for the null hypothesis of a constant difference of smoothed
#' Huber estimators, a confidence interval, and the EL estimator.
#'
#' @inheritParams EL.statistic
#' @param Y A numeric vector of data values.
#' @param mu A number specifying the null hypothesis value for the difference.
#'   Default is \code{0}.
#' @param scaleX The scale estimate of sample \code{X}. Default is \code{1}.
#' @param scaleY The scale estimate of sample \code{Y}. Default is \code{1}.
#' @param VX The asymptotic variance of the initial (nonsmooth) Huber estimator
#'   for sample \code{X}. Default is \code{2.046}.
#' @param VY The asymptotic variance of the initial (nonsmooth) Huber estimator
#'   for sample \code{Y}. Default is \code{2.046}.
#' @param k Tuning parameter for the Huber estimator. Default is \code{1.35}.
#'
#' @details
#' A common choice for a robust scale estimate (parameters \code{scaleX} and
#' \code{scaleY}) is the median absolute deviation (MAD).
#'
#' @return A list of class \code{"htest"} with components:
#' \describe{
#'   \item{\code{estimate}}{The empirical likelihood estimate for the difference
#'     of two smoothed Huber estimators.}
#'   \item{\code{conf.int}}{A confidence interval for the difference of two
#'     smoothed Huber estimators.}
#'   \item{\code{p.value}}{The p-value for the test.}
#'   \item{\code{statistic}}{The value of the test statistic.}
#'   \item{\code{method}}{The character string
#'     \code{"Empirical likelihood smoothed Huber estimator difference test"}.}
#'   \item{\code{null.value}}{The hypothesised difference \code{mu} under \eqn{H_0}.}
#'   \item{\code{data.name}}{A character string giving the names of the data.}
#' }
#'
#' @references
#' Valeinis, J. and Cers, E. Extending the two-sample empirical likelihood.
#' Preprint: \url{http://home.lu.lv/~valeinis/lv/petnieciba/EL_TwoSample_2011.pdf}
#'
#' Hampel, F., Hennig, C. and Ronchetti, E. A. (2011). A smoothing principle for
#' the Huber and other location M-estimators. \emph{Computational Statistics &
#' Data Analysis}, \strong{55}(1), 324--337.
#'
#' @author E. Cers, J. Valeinis
#'
#' @seealso \code{\link{EL.means}}
#'
#' @examples
#' X <- rnorm(100)
#' Y <- rnorm(100)
#' t.test(X, Y)
#' EL.means(X, Y)
#' EL.Huber(X, Y)
#'
#' @keywords nonparametric smooth htest
#' @export
EL.Huber <- function(X,
                     Y,
                     mu = 0,
                     conf.level = 0.95,
                     scaleX = 1,
                     scaleY = 1,
                     VX = 2.046,
                     VY = 2.046,
                     k = 1.35)
{
    call <- match.call()
    mm <- function(globals)
    {
        globals$Hsigma1 <- scaleX
        globals$Hsigma2 <- scaleY
        globals$Hconst <- k
        globals$HV1 <- VX
        globals$HV2 <- VY
        d <- deltasolve(globals, X, Y, 0)
        names(d) <- "Huber estimator difference"
        ci <- ELconf(globals, X, Y, 0, conf.level, delta = d)
        ll <- ELLR(globals, X, Y, mu, 0)
        names(ll) <- "-2 log LR"
        attr(ci, "conf.level") <- conf.level
        res <- list(
            statistic   = ll,
            parameter   = c(df = 1),
            p.value     = 1 - pchisq(ll, 1),
            estimate    = d,
            conf.int    = ci,
            null.value  = c("Huber difference" = mu),
            alternative = "two.sided",
            method      = "Two-sample empirical likelihood Huber estimator difference test",
            data.name   = paste(deparse(call$X), "and", deparse(call$Y))
        )
        class(res) <- "htest"
        res
    }
    EL.local(mm, method = "huber")
}

#### Globals setup
set.globals <- function(globals) {
    globals$lambda1      <- 0
    globals$lambda2      <- 0
    globals$wval1        <- NULL
    globals$wval2        <- NULL
    globals$lambdarange1 <- 0
    globals$lambdarange2 <- 0
    globals$degenerate   <- FALSE
    globals$precision    <- 1e-4
    globals
}

set.bw <- function(globals, bw) {
    globals$bw.X  <- 0
    globals$bw.Y  <- 0
    globals$bw    <- bw
    globals$initbw <- bw.init
    globals
}

bw.init <- function(globals, X, Y) {
    if ((globals$bw.X != 0) && (globals$bw.Y != 0))
        return(globals)
    if (is.function(globals$bw)) {
        globals$bw.X <- globals$bw(X)
        globals$bw.Y <- globals$bw(Y)
    } else {
        globals$bw.X <- globals$bw[1]
        globals$bw.Y <- globals$bw[2]
    }
    globals
}

set.kernel <- function(globals) {
    globals$krnl  <- kernel.norm
    globals$dkrnl <- kernel.dnorm
    globals$ikrnl <- kernel.inorm
    globals
}

kernel.norm  <- function(t)
    pnorm(t)
kernel.dnorm <- function(t)
    dnorm(t)
kernel.inorm <- function(t)
    qnorm(t)

set.method <- function(globals, type) {
    valid <- c("pp", "qq", "roc", "mean", "fdiff", "qdiff", "huber")
    if (!(type %in% valid))
        stop(paste(
            "Unknown EL method:",
            type,
            "\nRecognized EL methods are:",
            paste(valid, collapse = ", ")
        ))
    for (nm in c(
        "w1",
        "w2",
        "alpha1",
        "alpha2",
        "trange",
        "deltarange",
        "thetarange",
        "deltarange.opt",
        "thetarange.opt",
        "use.smoothing"
    )) {
        globals[[nm]] <- get(paste(nm, type, sep = "."))
    }
    globals
}

### Mean
w1.mean      <- function(globals, X, theta, delta, t)
    X - theta
w2.mean      <- function(globals, X, theta, delta, t)
    X - theta + delta
alpha1.mean  <- function(globals, X, theta, delta, t)
    - 1
alpha2.mean  <- function(globals, X, theta, delta, t)
    - 1
trange.mean  <- function(globals, X, Y)
    stop("Mean differences are not dependent on the parameter t.")
deltarange.mean <- function(globals, X, Y, t) {
    tr  <- c(min(X), max(X))
    tol <- 0.01 * (tr[2] - tr[1])
    c(tr[1] - max(Y) + tol, tr[2] - min(Y) - tol)
}
deltarange.opt.mean <- deltarange.mean
thetarange.mean <- function(globals, X, Y, delta, t) {
    tr  <- c(min(X), max(X))
    dr  <- c(min(Y) + delta, max(Y) + delta)
    tol <- 0.01 * (tr[2] - tr[1])
    c(max(tr[1], dr[1]) + tol, min(tr[2], dr[2]) - tol)
}
thetarange.opt.mean <- thetarange.mean
use.smoothing.mean  <- FALSE

### Huber estimator
### Huber
w1.huber <- function(globals, X, theta, delta, t) {
    globals$Hsigma <- globals$Hsigma1
    globals$HV <- globals$HV1
    Hsmooth(globals, (X - theta) / globals$Hsigma)
}
w2.huber <- function(globals, X, theta, delta, t) {
    globals$Hsigma <- globals$Hsigma2
    globals$HV <- globals$HV2
    Hsmooth(globals, (X - theta + delta) / globals$Hsigma)
}
alpha1.huber <- function(globals, X, theta, delta, t) {
    globals$Hsigma <- globals$Hsigma1
    globals$HV <- globals$HV1
    - DHsmooth(globals, (X - theta) / globals$Hsigma) / globals$Hsigma
}
alpha2.huber <- function(globals, X, theta, delta, t) {
    globals$Hsigma <- globals$Hsigma2
    globals$HV <- globals$HV2
    - DHsmooth(globals, (X - theta + delta) / globals$Hsigma) / globals$Hsigma
}
trange.huber       <- function(globals, X, Y)
    stop("Huber estimator differences are not dependent on the parameter t.")
deltarange.huber   <- deltarange.mean
deltarange.opt.huber <- deltarange.huber
thetarange.huber   <- thetarange.mean
thetarange.opt.huber <- thetarange.huber
use.smoothing.huber <- FALSE

Hsmooth <- function(globals, x) {
    x.l <- (x - globals$Hconst) / globals$HV
    x.r <- (x + globals$Hconst) / globals$HV
    pnorm(x.l) * (globals$Hconst - x) +
        pnorm(x.r) * (globals$Hconst + x) -
        globals$Hconst + globals$HV * (dnorm(x.r) - dnorm(x.l))
}
DHsmooth <- function(globals, x) {
    pnorm((x + globals$Hconst) / globals$HV) -
        pnorm((x - globals$Hconst) / globals$HV)
}

### ROC
w1.roc <- function(globals, X, theta, delta, t)
    globals$krnl((theta - X) / globals$bw.X) - (1 - delta)
w2.roc <- function(globals, X, theta, delta, t)
    globals$krnl((theta - X) / globals$bw.Y) - (1 - t)
alpha1.roc <- function(globals, X, theta, delta, t)
    globals$dkrnl((theta - X) / globals$bw.X) / globals$bw.X
alpha2.roc <- function(globals, X, theta, delta, t)
    globals$dkrnl((theta - X) / globals$bw.Y) / globals$bw.Y
trange.roc <- function(globals, X, Y)
    c(0.01, 0.99)
deltarange.roc <- function(globals, X, Y, t) {
    dY    <- globals$bw.Y * globals$ikrnl(t)
    diffs <- c(max(X) - min(Y) + dY, min(X) - max(Y) + dY)
    r     <- globals$krnl(diffs / globals$bw.X)
    c(min(r) + globals$precision, max(r) - globals$precision)
}
deltarange.opt.roc <- function(globals, X, Y, t) {
    kk <- 1 - ecdf(X)(c(
        quantile(Y, 1 - t) - 3 * globals$bw.X,
        quantile(Y, 1 - t) + 3 * globals$bw.X
    ))
    c(max(0.001, kk[2] - 0.001), min(0.999, kk[1] + 0.001))
}
thetarange.roc <- function(globals, X, Y, delta, t) {
    dX <- globals$bw.X * globals$ikrnl(delta)
    dY <- globals$bw.Y * globals$ikrnl(t)
    dd <- c(max(min(X) - dX, min(Y) - dY), min(max(X) - dX, max(Y) - dY))
    c(dd[1] + 0.01 * (dd[2] - dd[1]), dd[2] - 0.01 * (dd[2] - dd[1]))
}
thetarange.opt.roc <- function(globals, X, Y, delta, t) {
    extreme <- thetarange.roc(globals, X, Y, delta, t)
    tol     <- 2 * globals$bw.Y
    qY      <- quantile(Y, 1 - t)
    c(max(extreme[1], qY - tol), min(extreme[2], qY + tol))
}
use.smoothing.roc <- TRUE


### P-P plots
w1.pp <- function(globals, X, theta, delta, t)
    globals$krnl((theta - X) / globals$bw.X) - delta
w2.pp <- function(globals, X, theta, delta, t)
    globals$krnl((theta - X) / globals$bw.Y) - t
alpha1.pp <- function(globals, X, theta, delta, t)
    globals$dkrnl((theta - X) / globals$bw.X) / globals$bw.X
alpha2.pp <- function(globals, X, theta, delta, t)
    globals$dkrnl((theta - X) / globals$bw.Y) / globals$bw.Y
trange.pp <- function(globals, X, Y)
    c(0.01, 0.99)
deltarange.pp <- function(globals, X, Y, t)
    c(0.001, 0.999)
deltarange.opt.pp <- function(globals, X, Y, t) {
    kk <- ecdf(X)(c(
        quantile(Y, t) - 3 * globals$bw.X,
        quantile(Y, t) + 3 * globals$bw.X
    ))
    c(max(0.001, kk[1] - 0.001), min(0.999, kk[2] + 0.001))
}
thetarange.pp <- function(globals, X, Y, delta, t) {
    dX <- globals$bw.X * globals$ikrnl(delta)
    dY <- globals$bw.Y * globals$ikrnl(t)
    dd <- c(max(min(X) + dX, min(Y) + dY), min(max(X) + dX, max(Y) + dY))
    c(dd[1] + 0.01 * (dd[2] - dd[1]), dd[2] - 0.01 * (dd[2] - dd[1]))
}
thetarange.opt.pp <- function(globals, X, Y, delta, t) {
    dX      <- globals$bw.X * globals$ikrnl(delta)
    dY      <- globals$bw.Y * globals$ikrnl(t)
    dd      <- c(max(min(X) + dX, min(Y) + dY), min(max(X) + dX, max(Y) + dY))
    extreme <- c(dd[1] + 0.01 * (dd[2] - dd[1]), dd[2] - 0.01 * (dd[2] - dd[1]))
    tol <- 2 * globals$bw.Y
    qY <- quantile(Y, t)
    c(max(extreme[1], qY - tol), min(extreme[2], qY + tol))
}
use.smoothing.pp <- TRUE

### Q-Q plots
w1.qq <- function(globals, X, theta, delta, t)
    globals$krnl((delta - X) / globals$bw.X) - theta
w2.qq <- function(globals, X, theta, delta, t)
    globals$krnl((t - X) / globals$bw.Y) - theta
alpha1.qq <- function(globals, X, theta, delta, t)
    - 1
alpha2.qq <- function(globals, X, theta, delta, t)
    - 1
trange.qq <- function(globals, X, Y)
    c(min(Y), max(Y))
deltarange.qq <- function(globals, X, Y, t)
    c(min(X) + globals$precision, max(X) - globals$precision)
deltarange.opt.qq <- function(globals, X, Y, t) {
    extreme <- deltarange.qq(globals, X, Y, t)
    tol  <- globals$bw.Y
    tr   <- ecdf(Y)(c(t - tol, t + tol))
    upper <- quantile(X, tr[2]) + globals$bw.X
    lower <- quantile(X, tr[1]) - globals$bw.X
    c(max(lower, extreme[1]), min(upper, extreme[2]))
}
thetarange.qq <- function(globals, X, Y, delta, t) {
    ww <- w2.qq(globals, Y, 0, 0, t)
    tr <- c(min(ww), max(ww))
    ww <- w1.qq(globals, X, 0, delta, 0)
    dr <- c(min(ww), max(ww))
    c(max(tr[1], dr[1]) + globals$precision,
      min(tr[2], dr[2]) - globals$precision)
}
thetarange.opt.qq  <- thetarange.qq
use.smoothing.qq   <- TRUE

### F difference
w1.fdiff <- function(globals, X, theta, delta, t)
    globals$krnl((t - X) / globals$bw.X) - theta
w2.fdiff <- function(globals, X, theta, delta, t)
    globals$krnl((t - X) / globals$bw.Y) - theta - delta
alpha1.fdiff <- function(globals, X, theta, delta, t)
    - 1
alpha2.fdiff <- function(globals, X, theta, delta, t)
    - 1
trange.fdiff <- function(globals, X, Y)
    c(max(c(min(X), min(Y))), min(c(max(X), max(Y))))
deltarange.fdiff.prim <- function(globals, X, Y, t, prec) {
    ww <- w1.fdiff(globals, X, 0, 0, t)
    tr <- c(min(ww), max(ww))
    ww <- w2.fdiff(globals, Y, 0, 0, t)
    c(min(ww) - tr[2] + prec, max(ww) - tr[1] - prec)
}
deltarange.fdiff <- function(globals, X, Y, t)
    deltarange.fdiff.prim(globals, X, Y, t, globals$precision)
deltarange.opt.fdiff <- function(globals, X, Y, t)
    deltarange.fdiff.prim(globals, X, Y, t, 1e-7)
thetarange.fdiff <- function(globals, X, Y, delta, t) {
    ww  <- w1.fdiff(globals, X, 0, 0, t)
    tr  <- c(min(ww), max(ww))
    ww  <- w2.fdiff(globals, Y, 0, delta, t)
    tr2 <- c(min(ww), max(ww))
    c(max(tr2[1], tr[1]) + 1e-7, min(tr2[2], tr[2]) - 1e-7)
}
thetarange.opt.fdiff <- thetarange.fdiff
use.smoothing.fdiff  <- TRUE

### Quantile difference
w1.qdiff <- function(globals, X, theta, delta, t)
    globals$krnl((theta - X) / globals$bw.X) - t
w2.qdiff <- function(globals, X, theta, delta, t)
    globals$krnl((theta + delta - X) / globals$bw.Y) - t
alpha1.qdiff <- function(globals, X, theta, delta, t)
    globals$dkrnl((theta - X) / globals$bw.X) / globals$bw.X
alpha2.qdiff <- function(globals, X, theta, delta, t)
    globals$dkrnl((theta + delta - X) / globals$bw.Y) / globals$bw.Y
trange.qdiff <- function(globals, X, Y) {
    minlen <- min(length(X), length(Y))
    c(1 / minlen, 1 - 1 / minlen)
}
deltarange.qdiff <- function(globals, X, Y, t)
    c(min(Y) - max(X) + globals$precision,
      max(Y) - min(X) - globals$precision)
deltarange.opt.qdiff <- function(globals, X, Y, t) {
    extreme <- deltarange.qdiff(globals, X, Y, t)
    tol   <- 2 * globals$bw.Y
    qdiff <- quantile(Y, t) - quantile(X, t)
    c(max(extreme[1], qdiff - tol), min(extreme[2], qdiff + tol))
}
thetarange.qdiff <- function(globals, X, Y, delta, t) {
    dX <- globals$bw.X * globals$ikrnl(1 - t)
    dY <- globals$bw.Y * globals$ikrnl(1 - t)
    c(
        max(min(X) - dX, min(Y) - delta - dY) + globals$precision,
        min(max(X) - dX, max(Y) - delta - dY) - globals$precision
    )
}
thetarange.opt.qdiff <- function(globals, X, Y, delta, t) {
    extreme <- thetarange.qdiff(globals, X, Y, delta, t)
    tol <- 2 * globals$bw.X
    qX <- quantile(X, t)
    c(max(extreme[1], qX - tol), min(extreme[2], qX + tol))
}
use.smoothing.qdiff <- TRUE

EL.local <- function(fun, ..., bw = bw.nrd0, method) {
    globals <- list()
    globals <- set.method(globals, method)
    globals <- set.globals(globals)
    globals <- set.kernel(globals)
    globals <- set.bw(globals, bw)
    fun(globals, ...)
}

#### Calculations

ELconf <- function(globals,
                   X,
                   Y,
                   t,
                   p.level = 0.95,
                   delta = NULL) {
    globals <- globals$initbw(globals, X, Y)
    if (is.null(delta))
        delta <- deltasolve(globals, X, Y, t)
    if (is.na(delta))
        return(c(NA_real_, NA_real_))

    critval <- qchisq(p.level, 1)
    ELlim   <- function(delt)
        ELLR(globals, X, Y, delt, t) - critval

    drange <- globals$deltarange(globals, X, Y, t)
    if (drange[2] < drange[1])
        return(c(NA_real_, NA_real_))

    if (ELlim(delta) < 0) {
        lo <- if (ELlim(drange[1]) > 0) {
            uniroot(ELlim, c(drange[1], delta))$root
        } else {
            warning(paste("Could not find lower confidence bound for t =", t, "."),
                    call. = FALSE)
            drange[1]
        }
        # Upper bound
        hi <- if (ELlim(drange[2]) > 0) {
            uniroot(ELlim, c(delta, drange[2]))$root
        } else {
            warning(paste("Could not find upper confidence bound for t =", t, "."),
                    call. = FALSE)
            drange[2]
        }
    } else {
        warning(paste("Could not find estimate at t =", t, "."), call. = FALSE)
        lo <- hi <- delta
    }
    if (lo > hi) {
        warning(paste("Inverted confidence interval at t =", t, "; returning NA."),
                call. = FALSE)
        return(c(NA_real_, NA_real_))
    }
    c(lo, hi)
}

deltasolve <- function(globals, X, Y, t) {
    globals <- globals$initbw(globals, X, Y)
    dr <- globals$deltarange.opt(globals, X, Y, t)
    if (dr[2] < dr[1])
        return(NA_real_)
    if (dr[2] - dr[1] < 1e-7)
        return(dr[1])
    globals$tr         <- globals$thetarange
    globals$thetarange <- globals$thetarange.opt
    oo <- suppressWarnings(optimize(function(d)
        ELLR(globals, X, Y, d, t), dr))
    globals$thetarange <- globals$tr
    if (is.finite(oo$objective))
        oo$minimum
    else
        NA_real_
}

ELLR <- function(globals, X, Y, delta, t) {
    globals <- globals$initbw(globals, X, Y)
    ts      <- thetasolve(globals, X, Y, delta, t)
    globals <- ts$globals
    if (!globals$degenerate) {
        lr <- 2 * (sum(log(
            1 + globals$lambda1 * globals$wval1
        )) +
            sum(log(
                1 + globals$lambda2 * globals$wval2
            )))
        if (lr < 1e-7)
            0
        else
            lr
    } else {
        Inf
    }
}

thetasolve <- function(globals, X, Y, delta, t) {
    globals <- globals$initbw(globals, X, Y)
    len1 <- length(X)
    len2 <- length(Y)
    tgrad <- function(theta) {
        globals$degenerate <<- FALSE
        globals$wval1 <<- globals$w1(globals, X, theta, delta, t)
        globals$wval2 <<- globals$w2(globals, Y, theta, delta, t)
        alphaval1 <- globals$alpha1(globals, X, theta, delta, t)
        if (is.na(alphaval1[2]))
            alphaval1 <- rep.int(alphaval1, len1)
        alphaval2 <- globals$alpha2(globals, Y, theta, delta, t)
        if (is.na(alphaval2[2]))
            alphaval2 <- rep.int(alphaval2, len2)
        #' @useDynLib EL, .registration = TRUE
        res <- .C(
            "theta_equation",
            as.integer(len1),
            as.double(globals$wval1),
            as.double(alphaval1),
            as.integer(len2),
            as.double(globals$wval2),
            as.double(alphaval2),
            l1 = double(1),
            l2 = double(1),
            res = double(1)
        )
        globals$lambda1 <<- res$l1
        globals$lambda2 <<- res$l2
        res$res
    }
    tr <- globals$thetarange(globals, X, Y, delta, t)
    if (tr[1] > tr[2]) {
        globals$degenerate <- TRUE
        return(list(theta = NA, globals = globals))
    }
    theta <- suppressWarnings(try(uniroot(tgrad, tr)$root, silent = TRUE))
    list(theta = theta, globals = globals)
}

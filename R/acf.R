"acf" <-
function (x, lag.max = NULL, plot = FALSE, type = c("correlation", 
    "covariance", "partial")) 
{
    type <- match.arg(type)
    series <- deparse(substitute(x))
    x.freq <- frequency(as.ts(x))
    x <- as.matrix(x)
    sampleT <- nrow(x)
    #  or ? lag.max <- floor(10 * log10(sampleT))
    if (is.null(lag.max)) 
        lag.max <- floor(10 * (log10(sampleT) - log10(ncol(x))))
    lag.max <- min(lag.max, sampleT - 1)
    xb <- sweep(x, 2, apply(x, 2, mean))
    lag <- matrix(1, ncol(x), ncol(x))
    lag[lower.tri(lag)] <- -1
    if (type == "partial") {
        acf <- ar(x, order.max = lag.max)$partialacf
        lag <- outer(1:lag.max, lag/x.freq)
    }
    else {
        acf <- array(NA, c(lag.max + 1, ncol(x), ncol(x)))
        for (i in 0:lag.max) {
            Om <- (t(xb[(i + 1):sampleT, , drop = FALSE]) %*% 
                xb[1:(sampleT - i), , drop = FALSE])/sampleT
            if (type == "correlation") {
                # nrow above for univariate case
                if (i == 0) 
                  Om0 <- diag(1/sqrt(diag(Om)), nrow = nrow(Om))
                Om <- Om0 %*% Om %*% Om0
            }
            acf[i + 1, , ] <- Om
        }
        lag <- outer(0:lag.max, lag/x.freq)
    }
    acf.out <- structure(.Data = list(acf = acf, type = type, 
        n.used = sampleT, lag = lag, series = series, snames = colnames(x)), 
        class = "acf")
    if (plot) 
        plot.acf(acf.out)
    else return(acf.out)
}
"plot.acf" <-
function (x, ci = 0.95, type = "h", xlab = "Lag", ylab = NULL, 
    ylim = NULL, main = NULL, ...) 
{
    opar <- NULL
    on.exit(par(opar))
    nser <- ncol(x$lag)
    opar <- c(opar, par(mfrow = rep(min(nser, 5), 2)))
    if (is.null(ylab)) 
        ylab <- switch(x$type, correlation = "ACF", covariance = "ACF", 
            partial = "Partial ACF")
    if (is.null(snames <- x$snames)) {
        snames <- if (nser == 1) 
            paste("Series", x$series)
        else paste("Series", 1:nser)
    }
    with.ci <- (ci != 0) && (x$type != "covariance")
    for (i in 1:nser) for (j in 1:nser) {
        clim <- if (with.ci) 
            qnorm((1 + ci)/2)/sqrt(x$n.used)
        else c(0, 0)
        if (is.null(ylim)) {
            ymin <- min(c(x$acf[, i, j], -clim))
            ymax <- max(c(x$acf[, i, j], clim))
            ylim <- c(ymin, ymax)
        }
        plot(x$lag[, i, j], x$acf[, i, j], type = type, xlab = xlab, 
            ylab = ylab, ylim = ylim, ...)
        if (with.ci) 
            abline(h = c(clim, -clim), col = "red", lty = 2)
        abline(h = 0)
        if (!is.null(main)) 
            title(main)
        else if (i == j) 
            title(snames[i])
        else title(paste(snames[i], "&", snames[j]))
    }
}

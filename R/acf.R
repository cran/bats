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
        n.used = sampleT, lag = lag, series = series), class = "acf")
    if (plot) 
        plot.acf(acf.out)
    else return(acf.out)
}
"plot.acf" <-
function (x, coverage = 0.95, ...) 
{
    opar <- NULL
    on.exit(par(opar))
    nser <- dim(x$lag)[2]
    m <- min(nser, 5)
    opar <- c(opar, par(mfrow = c(m, m)))
    ylab <- switch(x$type, correlation = "ACF", covariance = "ACF", 
        partial = "Partial ACF")
    if (nser > 1) 
        if (is.null(dimnames(x)) || is.null(dimnames(x)[[2]])) 
            names <- paste("Series", 1:nser)
        else names <- dimnames(x)[[2]]
    for (i in 1:nser) for (j in 1:nser) {
        if (x$type == "correlation" || x$type == "partial") {
            clim <- qnorm((1 + coverage)/2)/sqrt(x$n.used)
            ymin <- min(c(x$acf[, i, j], -clim))
            ymax <- max(c(x$acf[, i, j], clim))
            plot(x$lag[, i, j], x$acf[, i, j], type = "h", xlab = "Lag", 
                ylab = ylab, ylim = c(ymin, ymax), ...)
            abline(h = c(clim, -clim), col = "red", lty = 2)
        }
        else {
            plot(x$lag[, i, j], x$acf[, i, j], type = "h", xlab = "Lag", 
                ylab = ylab, ...)
        }
        abline(h = 0)
        if (nser == 1) 
            title(paste("Series:", x$series))
        else {
            if (i == j) 
                title(names[i])
            else title(paste(names[i], "&", names[j]))
        }
    }
}

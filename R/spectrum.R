"spectrum" <- 
function (x, method = c("pgram", "ar"), plot = TRUE, ...)
{
    # Wrapper function for spec methods
    # (need to implement spec.ar)
    method <- match.arg(method)
    switch(method,
        pgram = spec.pgram(x, plot = plot, ...),
        ar = stop("spec.ar not implemented yet"))
}
"spec.pgram" <-
function (x, spans = 1, taper = 0.1, demean = FALSE, detrend = TRUE, 
    pad = 0, plot = FALSE) 
{
    # Estimate spectral density from smoothed periodogram.
    # (we could farm out padding and tapering to other functions)
    #
    series <- deparse(substitute(x))
    xfreq <- frequency(x)
    x <- as.matrix(x)
    N <- nrow(x)
    nser <- ncol(x)
    if (detrend) {
        t <- 1:N - (N + 1)/2
        sumt2 <- N * (N^2 - 1)/12
        for (i in 1:ncol(x)) {
            x[, i] <- x[, i] - mean(x[, i]) - sum(x[, i] * t) * t/sumt2
	}
        
    }
    else if (demean) {
        x <- sweep(x, 2, apply(x, 2, mean))
    }
    if (taper > 0.5 || taper < 0) 
        stop("taper must be between 0 and 0.5")
    else if (taper > 0) {
        w <- rep(1, N)
        n <- max(round(N * taper), 1)
        w[N:(N - n + 1)]  <- w[1:n] <- sin(((1:n - 0.5) * pi)/(2 * n))^2
        x <- x * w
    }
    if (pad > 0) {
        x <- rbind(x, matrix(0, nrow = N * pad, ncol = ncol(x)))
        N <- nrow(x)
    }
    NewN <- nextn(N)
    x <- rbind(x, matrix(0, nrow = (NewN - N), ncol = ncol(x)))
    N <- nrow(x)
    Nspec <- floor(N/2)
    freq <- seq(from = xfreq/N, by = xfreq/N, length = Nspec)
    xfft <- mvfft(x)[2:(N - 1), , drop = FALSE]
    pgram <- array(NA, dim = c(nrow(x) - 2, ncol(x), ncol(x)))
    for (i in 1:ncol(x)) {
        for (j in 1:ncol(x)) {
            pgram[, i, j] <- xfft[, i] * Conj(xfft[, j])/(N*xfreq)
        }
    }
    filter.list <- vector("list", length(spans))
    for (i in 1:length(spans)) {
        m <- floor(spans[i]/2)
        spans[i] <- 2 * m + 1
        filter.list[[i]] <- if (m > 0) 
            c(0.5, rep(1, 2 * m - 1), 0.5)/(2 * m)
        else 1
    }
    filter <- filter.list[[1]]
    if (length(spans) > 1) 
        for (i in 2:length(spans)) filter <- convolve(filter.list[[i]], 
            filter, type="open")
    if (length(filter) > 1) {
        ndiff <- nrow(pgram) - length(filter)
        m <- floor(length(filter)/2)
        if (ndiff < 0) 
            stop("filter too long!")
        else for (i in 1:ncol(x)) for (j in 1:ncol(x)) {
            pgram[, i, j] <- convolve(pgram[, i, j], c(filter[(m + 
                1):(2 * m + 1)], rep(0, ndiff), filter[1:m]))
        }
    }
    df <- 2/sum(filter^2)
    m <- floor(length(filter)/2)
    bandwidth <- sqrt(sum((1/12 + (0:(2 * m) - m)^2) * filter)) * xfreq/nrow(x)
    spec <- matrix(NA, nrow = Nspec, ncol = nser)
    for (i in 1:nser) spec[, i] <- Re(pgram[1:Nspec, i, i])
    if (nser == 1) {
        coh <- phase <- NULL
        spec <- drop(spec)
    }
    else {
        coh <- phase <- matrix(NA, nrow = Nspec, ncol = nser * 
            (nser - 1)/2)
        for (i in 1:(nser - 1)) {
            for (j in (i + 1):nser) {
                coh[, i + (j - 1) * (j - 2)/2] <- Mod(pgram[1:Nspec, 
                  i, j])^2/(spec[, i] * spec[j])
                ph <- Arg(pgram[1:Nspec, i, j])
                dph <- diff(ph)
                ph <- c(0, sign(dph) * (abs(dph)%%pi)) + ph[1]
                phase[, i + (j - 1) * (j - 2)/2] <- ph
            }
        }
    }
    spec <- 10 * log10(spec)
    spg.out <- list(freq = freq, spec = spec, coh = coh, phase = phase, 
        spans = spans, filter = filter, df = df, bandwidth = bandwidth, 
        n.used = nrow(x), series = series, method = ifelse(length(filter) > 
            1, "Smoothed Periodogram", "Raw Periodogram"), taper = taper, 
        pad = pad, detrend = detrend, demean = demean)
    class(spg.out) <- "spec"
    if(plot)
	plot.spec(spg.out)
    else return(spg.out)
}
"plot.spec" <-
function (x, add = FALSE, ci = 0.95, xlab = "frequency", 
    ylab = "spectrum (dB)", type = "l", main = NULL, sub = NULL, ...) 
{
    matplot(x$freq, x$spec, xlab = xlab, ylab = ylab, type = type, 
        add = add, ...)
    if (ci == 0 || add) {
        #No confidence limits
        ci.text <- ""
    }
    else {
        # The position of the error bar has no meaning: only the width
        # and height. It is positioned in the top right hand corner.
        #
        conf.lim <- spec.ci(x, coverage = ci)
        conf.y <- max(x$spec) - conf.lim[2]
        conf.x <- max(x$freq) - x$bandwidth
        lines(rep(conf.x, 2), conf.y + conf.lim)
        lines(conf.x + c(-0.5, 0.5) * x$bandwidth, rep(conf.y, 
            2))
        ci.text <- paste("95% C.I. is (", paste(format(conf.lim, 
            digits = 3), collapse = ","), ")dB")
    }
    if (is.null(main)) 
        main <- paste(paste("Series:", x$series), x$method, sep = "\n")
    if (is.null(sub)) 
        sub <- paste("bandwidth=", format(x$bandwidth, digits = 3), 
            ci.text)
    title(main = main, sub = sub)
    invisible(x)
}
"spec.ci" <-
function (spec.obj, coverage = 0.95) 
{
    # A utility function for plot.spec which calculates the confidence
    # interval (centred around zero). We use a conditional argument to
    # ensure that the ci always contains zero.
    #
    if (coverage < 0 || coverage >= 1)
        stop("coverage probability out of range [0,1)")
    df <- spec.obj$df
    limits <- numeric(2)
    upper.quantile <- 1 - (1 - coverage) * (1 - pchisq(df, df))
    lower.quantile <- (1 - coverage) * pchisq(df, df)
    -10 * log10(qchisq(c(upper.quantile, lower.quantile), df)/df)
}

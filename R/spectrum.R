"plot.spec" <-
function (x, ...) 
{
    matplot(x$freq, x$spec, xlab = "frequency", ylab = "spectrum", 
        type = "l", xlim = c(0, 0.5), ...)
    title(main = paste("Series:", x$series), sub = paste("bandwidth=", 
        format(x$bandwidth, digits=4)))
    mtext(x$method, 3, 0.5)
}
"spec.pgram" <-
function (x, spans = 1, taper = 0.1, demean = FALSE, detrend = TRUE, 
    pad = 0, plot = FALSE) 
{
    series <- deparse(substitute(x))
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
        w[1:n] <- sin(((1:n - 0.5) * pi)/(2 * n))^2
        w[N:(N - n + 1)] <- w[1:n]
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
    freq <- seq(from = 1/N, by = 1/N, length = Nspec)
    xfft <- mvfft(x)[2:(N - 1), , drop = FALSE]
    pgram <- array(NA, dim = c(nrow(x) - 2, ncol(x), ncol(x)))
    for (i in 1:ncol(x)) {
        for (j in 1:ncol(x)) {
            pgram[, i, j] <- xfft[, i] * Conj(xfft[, j])/N
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
        for (i in 2:length(spans)) filter <- convolve.open(filter.list[[i]], 
            filter)
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
    bandwidth <- sqrt(sum((1/12 + (0:(2 * m) - m)^2) * filter))/nrow(x)
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
    class(spg.out) <- "spec.pgram"
    if(plot)
	plot.spec(spg.out)
    else return(spg.out)
}

"convolve" <-
function (x, y, conj = TRUE) 
{
    n <- length(x)
    if (length(y) != n) 
        stop("length mismatch in convolution")
    z <- fft(fft(x) * (if (conj) 
        Conj(fft(y))
    else fft(y)), inv = TRUE)/n
    if (is.real(x) && is.real(y)) 
        Re(z)
    else z
}
"convolve.open" <-
function (x, y, conj = TRUE) 
{
    nx <- length(x)
    ny <- length(y)
    n <- nx + ny - 1
    x <- c(rep(0, ny - 1), x)
    y <- c(y, rep(0, nx - 1))
    z <- fft(fft(x) * (if (conj) 
        Conj(fft(y))
    else fft(y)), inv = TRUE)/n
    if (is.real(x) && is.real(y)) 
        Re(z)
    else z
}

\name{spec.pgram}
\alias{spec.pgram}
\title{Estimate spectral density of a time series from smoothed periodogram}
\usage{
spec.pgram(x, spans=1, taper=0.1, demean=FALSE, detrend=TRUE,
	pad=0, plot=FALSE)
}
\arguments{
 \item{x}{A time series.}
 \item{spans}{Vector of odd integers giving the widths of modified Daniell
  smoothers to be used to smooth the periodogram.}
 \item{taper}{Proportion of data to taper.  A split cosine bell taper is
  applied to this proportion of the data at the beginning and end of
  the series.}
 \item{demean}{logical. If \code{TRUE}, subtract the mean of the series.}
 \item{detrend}{logical. If \code{TRUE}, remove a linear trend from the
    series. This will also remove the mean.}
 \item{pad}{Proportion of data to pad. Zeros are added to the end of
  the series to increase its length by the proportion \code{pad}.}
 \item{plot}{Plot the periodogram.}
} \description{
\code{spec.pgram} calculates the periodogram using a fast Fourier
transform, and optionally smooths the result with a series of modified
Daniell smoothers (moving averages giving half weight to the end values).
}
\details{
The raw periodogram is not a consistent estimator of the spectral density,
but adjacent values are asymptotically independent. Hence a consistent
estimator can be derived by smoothing the raw periodogram, assuming that
the spectral density is smooth.
}
\value{
   An object of class \code{\link{spec}} with the following additional
   list elements.
   \item{spans}{The \code{spans} argument. Values which are not odd
      integers are converted to the next highest odd integer.}
   \item{filter}{The convolution of the modified Daniell smoothers
      which is used to smooth the raw periodogram.}
   \item{df}{The distribution of the spectral density estimate can
      be approximated by a chi square distribution with \code{df}
      degrees of freedom.}
   \item{bandwidth}{}
   \item{taper}{The value of the \code{taper} argument}
   \item{pad}{The value of the \code{pad} argument}
   \item{detrend}{The value of the \code{detrend} argument}
   \item{demean}{The value of the \code{demean} argument}
}
\note{The series will be automatically padded with zeros until the
series length is a highly composite number in order to help the Fast
Fourier Transform. This is done automatically and is not controlled by
the \code{pad} argument.}

%\references{ ~put references to the literature/web site here ~ }
\author{Martyn Plummer}
\seealso{\code{\link{spectrum}},\code{\link{plot.spec}}, \code{\link{fft}}}
%\examples{}
\keyword{ts}%-- one or more ...

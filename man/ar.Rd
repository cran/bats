\name{ar}
\alias{ar}
\title{Fit autoregressive models to time series}
\usage{
ar(x, aic=TRUE, order.max=NULL, method=c("yule-walker", "burg"))
}
\arguments{
 \item{x}{A time series}
 \item{aic}{Logical flag. If \code{TRUE} then the Akaike Information
     Criterion is used to choose the order of the autoregressive model.}
 \item{order.max}{Maximum order of model to fit. Defaults to 10*log10(N)
     where N is the number of observations.}
 \item{method}{Character string giving the method used to fit the model.
     Must be one of \code{yule-walker} or \code{burg} (the first few
     characters are sufficient). Defaults to \code{yule-walker}.}
}
\description{
This is just a wrapper for the functions \code{\link{ar.yw}} and
\code{\link{ar.burg}}.
}
\value{
  A list with the following elements:
  \item{order}{The order of the fitted model.  This is chosen by
     minimizing the AIC if \code{aic=TRUE}, otherise it is \code{order.max}.}
  \item{ar}{Estimated autoregression coefficients for the fitted model.
     A three dimensional array.}
  \item{var.pred}{The prediction variance: the portion of the
     variance of the time series that is not explained by the
     autoregressive model.}
  \item{aic}{The value of the \code{aic} argument.}
  \item{n.used}{The number of observations in the time series.}
  \item{order.max}{The value of the \code{order.max} argument.}
  \item{partialacf}{The estimate of the partial autocorrelation function
     up to lag \code{order.max}}
  \item{resid}{residuals from the fitted model, conditioning on the
     first \code{order} observations. The first
     \code{order} residuals are set to \code{NA}. 
  \item{method}{The value of the \code{method} argument.}
  \item{series}{The name of the time series}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Martyn Plummer}
%\note{ ~~further notes~~ }
\seealso{\code{\link{ar.yw}}}

%\examples{}
\keyword{ts}

\name{ar.yw}
\alias{ar.yw}
\title{Fit autoregressive model by solving Yule-Walker equations}
\usage{
ar.yw(x, aic=TRUE, order.max=NULL)
}
\arguments{
 \item{x}{A time series.}
 \item{aic}{Logical flag. If \code{TRUE} then use Akaike Information
    criterion to select order of model.}
 \item{order.max}{The maximum order of autoregressive model.}
}
\description{

}
\value{
   See \code{\link{ar}}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Martyn Plummer}
\note{\code{ar.yw} is called by \code{acf} to calculate the
partial autocorrelation function.}

\seealso{ \code{\link{ar}}, \code{\link{ar.burg}}, \code{\link{acf}} }
%\examples{}
\keyword{ts}

\name{summary.ICCfit}
\alias{summary.ICCfit}
\alias{print.summary.ICCfit}
\title{
Summary for ICCfit Objects
}
\description{
The function summarizes the results of ICC estimation.
}
\usage{
\method{summary}{ICCfit}(object, ...)
}
\arguments{
  \item{object }{
object of the class ICCfit (output of the \code{\link{intervalICC}} function)
}
\item{... }{additional arguments passed to the function (they do not affect the summary produced)
}
}
\details{
For more details about ICC estimation and output values shortly described below, please refer to the documentation for \code{\link{intervalICC}}.
}
\value{
An object of class "summary.ICCfit". The object is a list with the components:
  \item{estimates }{a data frame containing maximum likelihood estimates for ICC, mean and variance components}
  \item{loglikelihood }{log-likelihood evaluated at maximum likelihood estimates}
}
\references{
Kovacic J, Varnai VM. Reproducibility of grouped data assessed by intraclass correlation coefficient. Accepted for publication in \emph{Epidemiology}.
}
\author{
Jelena Kovacic \email{jkovacic@imi.hr}
}
\seealso{
\code{\link{intervalICC}}
}
\examples{
# Example with 6 predefined classes (grouped data)
classes <- 1:6
class.limits <- cbind(classes-0.5,classes+0.5)
r1 <- sample(classes,30,replace=TRUE)
r2 <- sample(classes,30,replace=TRUE)
icc.est <- intervalICC(r1,r2,predefined.classes=TRUE,classes,class.limits)
summary(icc.est)
}

\name{iRepro-package}
\alias{iRepro-package}
\alias{iRepro}
\docType{package}
\title{
Reproducibility for Interval-Censored Data
}
\description{
Calculates intraclass correlation coefficient (ICC) for assessing reproducibility of interval-censored data with two repeated measurements (Kovacic and Varnai (2014) <doi:10.1097/EDE.0000000000000139>). ICC is estimated by maximum likelihood from model with one fixed and one random effect (both intercepts). Help in model checking (normality of subjects' means and residuals) is provided.
}
\details{
\tabular{ll}{
Package: \tab iRepro\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2018-06-12\cr
License: \tab GPL-3\cr
}
}
\author{
Jelena Kovacic

Maintainer: Jelena Kovacic <jkovacic@imi.hr>
}
\references{
Kovacic J, Varnai VM. Intraclass correlation coefficient for grouped data. \emph{Epidemiology} 2014;25(5):769--770.
}
\examples{
# Data generation (grouped data)
classes <- 1:6
class.limits <- cbind(classes-0.5,classes+0.5)
r1 <- sample(classes,100,replace=TRUE) # first measurement
r2 <- sample(classes,100,replace=TRUE) # second measurement

summary(intervalICC(r1,r2,predefined.classes=TRUE,classes,class.limits)) # ICC estimation
}

\name{PhenoEvoData_run4}
\alias{PhenoEvoData_run4}
\docType{data}
\title{
A small PhenoEvoDataset for example use - a single run, subsampled timepoints.
}
\description{
This dataset is a small sample dataset, a subset of the main PhenoEvoData example dataset available with our tutorial online. It contains data from only a single experiment (run4) and from only every 10th timepoint. It is useful for testing functions for plotting timecourses, such as switch.rate.timecourse() and barcode.timecourse().
PhenoEvoData_run4 was generated from PhenoEvoData using the following:
PhenoEvoData_run4<-subset(PhenoEvoData, run.number==4 & step %in% seq(0, 1000, 10))
PhenoEvoData_run4$step<-seq(1, length(PhenoEvoData_4$step))
}
\usage{data("PhenoEvoData_run4")}
\format{
  A data frame with 101 observations on the following 12 variables.
  \describe{
    \item{\code{run.number}}{a numeric vector}
    \item{\code{toxin.conc}}{a numeric vector}
    \item{\code{env.noise}}{a numeric vector}
    \item{\code{step}}{a numeric vector}
    \item{\code{count.turtles}}{a numeric vector}
    \item{\code{mean.toxin}}{a numeric vector}
    \item{\code{degrade.rate}}{a character vector}
    \item{\code{switch.rate}}{a character vector}
    \item{\code{response.error}}{a character vector}
    \item{\code{barcode}}{a character vector}
    \item{\code{generation}}{a character vector}
    \item{\code{xydr}}{a character vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(PhenoEvoData_run4)
switch.rate.timecourse(PhenoEvoData_run4, 4)
barcode.timecourse(PhenoEvoData_run4, 4)
}
\keyword{datasets}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addMarkerInformation.R
\name{addSNPdetails}
\alias{addSNPdetails}
\title{Add SNP details to the blast output}
\usage{
addSNPdetails(
  blast.path,
  reference.bases,
  probe.name,
  genome.name,
  output.path = "Mapping_results",
  pident.filter = 90,
  coverage.filter = 50
)
}
\arguments{
\item{reference.bases}{Path to the reference bases obtained from samtools}

\item{probe.name}{name of the probe file being analysed}

\item{genome.name}{name of the reference genome being compared}

\item{output.path}{path were the results are saved}

\item{coverage.filter}{double value to filter hits with small coverage}

\item{blast.results}{Dataframe or path to blast results returned by the `processBlastResults()` function}
}
\value{
Returns a dataframe with the updated blast results
}
\description{
This function is used to add more details to the blast output returned by the `processBlastResults()` function
}
\author{
David Chisanga
}
\keyword{allele}
\keyword{allele,}
\keyword{alternative}
\keyword{genome,}
\keyword{reference}

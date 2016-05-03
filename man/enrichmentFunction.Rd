\name{enrichmentFunction}
\alias{enrichmentFunction}

\title{
    Identification of statistically significant GO terms
}

\description{
    Using the hypergenomic distribution test and FDR correction procedure, this function can identify a set of
    statistically significant GO terms.
}

\arguments{
  \item{annRef}{
      \code{annRef} is the annotated reference genes and their annotated terms.
  }  
  \item{annInterest}{
       \code{annInterest} is the annotated interesting genes and their annotated terms.
  }
  \item{method}{
       \code{method} is the fdr method. GO-function provides three p value correction methods: "bonferroni", "BH" 
       and "BY". The default is "BY".
  }
  \item{fdrth}{
       \code{fdrth} is the fdr cutoff to identify statistically significant GO terms. The default is 0.05.
  }
}

\value{
     The function returns a list contained two values: The first value is a table contained the information of all
     statistically significant terms under a FDR control level and the second value is a table contained the information 
     of all annotated GO terms.
}

\author{
    Jing Wang
}

\seealso{
    \code{\link{GOFunction}}
}

\keyword{ methods }
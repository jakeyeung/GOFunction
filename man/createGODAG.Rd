\name{createGODAG}
\alias{createGODAG}


\title{
     Creation of GO DAG stucture for statistically significant GO terms
}

\description{
     To plot the relationship between statistically significant GO terms, this function creates a GO DAG structure 
     for these terms.
}

\arguments{
  \item{sigNodes}{
      \code{sigNodes} is the statistically significant GO terms found by "enrichmentFunction" function. 
  }
  \item{ontology}{
      The default \code{ontology} is "BP" (Biological Process). The "CC" (Cellular Component) and "MF" (Molecular 
      Function) ontologies can also be used.
  }
}

\value{
   This function returns a object of 'graphNEL' class.
}

\author{
   Jing Wang
}

\note{
   This function simulates the related program in TopGO (Alexa, A. et al. (2006) Improved scoring of functional groups
   from gene expression data by decorrelating GO graph structure. Bioinformatics, 22, 1600-1607). 
}

\seealso{
     \code{\link{GOFunction}}
     \code{\link{enrichmentFunction}}
}

\keyword{ methods }
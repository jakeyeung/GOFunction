\name{localRedundancy}
\alias{localRedundancy}

\title{
    Treating for local redundancy between GO terms with ancient-offspring relationship
}
\description{
    This function uses some explicit rules to treat local redundancy (see details).
}

\arguments{
  \item{sigTerm}{
       \code{sigTerm} is the statistically significant terms found by the enrichmentFunction function.
  }
  \item{generalAnn}{
       \code{generalAnn} is the mapping between gene and GO terms.
  }
  \item{sigTermRelation}{
        \code{sigTermRelation} is the GO DAG relationship between the statistically significant terms.
  }
  \item{annRef}{
        \code{annRef} is the annotated reference genes and their annotated terms.
  }
  \item{annInterest}{
        \code{annInterest} is the annotated interesting genes and their annotated terms.
  }
  \item{ppth}{
        \code{ppth} is the significant level to test whether the remaining genes of the ancestor term are enriched 
        with interesting genes after removing the genes in its significant offspring terms. The default is 0.05.
  }
  \item{pcth}{
       \code{pcth} is the significant level to test whether the frequency of interesting genes in the offspring 
       terms are significantly different from that in the ancestor term. The default is 0.05. 
  }
}


\details{
    This function uses three principles to treat local redundancy between GO terms with ancient-offspring relationship: 
    1. if the frequency of interesting genes in the remaining genes of the ancestor term is not significantly higher 
    than the random background, and at the same time, it is lower than the frequency of its significant offspring terms, 
    GO-function only remains the offspring terms; Otherwise, GO-function remains ancient term. In this situation, 
    2. if the offspring terms are no different with the ancient term, GO-function only remains the ancestor term; 
    otherwise, 3. GO-function remains both ancient term and its offspring terms.
}


\value{
    The function returns a set of remained terms without local redundancy.
}

\author{
   Jing Wang
}

\seealso{
    \code{\link{GOFunction}}
    \code{\link{enrichmentFunction}}
}

\keyword{ methods }
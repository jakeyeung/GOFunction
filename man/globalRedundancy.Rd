\name{globalRedundancy}
\alias{globalRedundancy}

\title{
     Treating for glocal redundancy between GO terms which share multiple function genes but have no ancient-offspring
     relationship 
}
\description{
     For a pair of terms with overlapping genes, this function can identify whether their significance might be simply
     introduced by the overlapping genes (see details).
}

\arguments{
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
  \item{sigTermRedun}{
     \code{sigTermRedun} is the remained statistically signifcant terms after treating local redundancy.
  }
  \item{poth}{
     \code{poth} is the significant level to test whether the overlapping genes of one term is significantly
     different from the non-overlapping genes of the term. The default is 0.05.
  }
  \item{peth}{
    \code{peth} is the significant level to test whether the non-overlapping genes of one term is enriched with
    interesting genes. The default is 0.05.
  }
}

\details{
     This function uses the following principles to treat for glocal redundancy between GO terms which share multiple 
     function genes but have no ancient-offspring relationship. For a pair of overlapping terms, GO-function remains 
     one term only if (1) there is additional evidence that the non-overlapping genes of a term may be related to the
     disease if the frequency of interesting genes in these genes is not significantly different from that in the
     overlapping genes or is significantly higher than the random background; and (2) there is no such evidence for 
     the non-overlapping genes of another term. In all the other situations, GO-function remains both terms.  
}
\value{
    The function returns a set of remained terms with evidence that their significance should not 
    be simply due to the overlapping genes.
}

\author{
     Jing Wang
}

\seealso{
     \code{\link{GOFunction}}
     \code{\link{localRedundancy}}
}

\keyword{ methods }
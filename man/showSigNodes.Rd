\name{showSigNodes}
\alias{showSigNodes}

\title{
    Drawing the GO DAG structure of statistically significant terms
}
\description{
    According to the GO DAG structure created by the "createGODAG" function, this function draws the corresponding plot.
}

\arguments{
  \item{DAG}{
     \code{DAG} is a object of 'graphNEL' class which contains the GO DAG structure of statistically significant GO terms.  
  }
  \item{sigTerm}{
     \code{sigTerm} is the GO terms found by enrichmentFunction function.
  }
  \item{sigTerm_Local}{
     \code{sigTerm_Local} is the remained statistically signifcant terms after treating local redundancy.
  }
  \item{sigTerm_Global}{
     \code{sigTerm_Glocal} is the remained statistically signifcant terms after treating global redundancy.
  }
  \item{dagTermInfo}{
     \code{dagTermInfo} is the information of all related GO terms in GO DAG.
  }
  \item{bmpSize}{
     \code{bmpSize} is the width and height of the plot of GO DAG for all statistically significant terms. GO-function
     set the default width and height of the plot as 2000 pixels in order to clearly show the GO DAG structure. If the 
     GO DAG is very complexity, the user should increase \code{bmpSize}. Note: If there is an error at the step of 
     "bmp(filename, width = 2000, ..." when running GO-function, the user should decrease \code{bmpSize}.
  }
  \item{filename}{
    \code{filename} is the name of the files saving the table and the GO DAG of all statistically significant terms.
  }

}

\value{
   This function saves the structure of GO DAG for all statistic significant terms into a plot (e.g. "sigTerm.bmp") in 
   the current working folder. In this plot, "circle" represents terms removed after local redundancy treatment, 
   "box" represents terms removed after treating global redundancy and "rectangle" represents the remained terms 
   with evidence that their significance should not be simply due to the overlapping genes. The dfferent 
   color shades represent the  adjusted p values of the terms.
}

\author{
    Jing Wang
}

\seealso{
     \code{\link{GOFunction}}
     \code{\link{enrichmentFunction}}
     \code{\link{localRedundancy}}
     \code{\link{globalRedundancy}}
     \code{\link{createGODAG}}
}

\keyword{ methods }
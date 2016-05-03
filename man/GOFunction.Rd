\name{GOFunction}
\alias{GOFunction}

\title{ main function of the GO-function package}

\description{
   The \code{GOFunction} function is the main function of the GO-function package and can generate a set of biologically
   relevant GO terms.
}
\usage{
GOFunction(interestGenes, refGenes, organism = "org.Hs.eg.db", ontology = "BP", fdrmethod = "BY", fdrth = 0.05, ppth =   0.05, pcth = 0.05, poth = 0.05, peth = 0.05, bmpSize = 2000, filename = "sigTerm")
}

\arguments{
  \item{interestGenes}{
    \code{interestGenes} is a set of interesting genes (e.g. differential expressed genes), which should be 
    denoted using the Entrez gene ID.
  }
  \item{refGenes}{
    \code{refGenes} is the background genes corresponding to the interesting genes, which should be denoted using 
    the Entrez gene ID.
  }
  \item{organism}{
     The GO-function package can be currently applied to analyse data for 18 organisms and the user should install the 
     corresponding gene annotation package when analysing data for these organisms. The 18 organisms and the 
     corresponding packages are as follows: Anopheles "org.Ag.eg.db", Bovine "org.Bt.eg.db", Canine "org.Cf.eg.db", 
     Chicken "org.Gg.eg.db", Chimp "org.Pt.eg.db", E coli strain K12 "org.EcK12.eg.db", E coli strain Sakai
     "org.EcSakai.eg.db", Fly "org.Dm.eg.db", Human "org.Hs.eg.db", Mouse "org.Mm.eg.db", Pig "org.Ss.eg.db", Rat
     "org.Rn.eg.db", Rhesus "org.Mmu.eg.db", Streptomyces coelicolor "org.Sco.eg.db", Worm "org.Ce.eg.db", Xenopus
     "org.Xl.eg.db", Yeast "org.Sc.sgd.db", Zebrafish "org.Dr.eg.db". The default \code{organism} is "org.Hs.eg.db" 
     (Human).
  }
  \item{ontology}{
       The default \code{ontology} is "BP" (Biological Process). The "CC" (Cellular Component) and "MF" (Molecular 
       Function) ontologies can also be used.
  }
  \item{fdrmethod}{
     GO-function provides three p value correction methods: "bonferroni", "BH" and "BY". The default \code{fdrmethod} 
     is "BY".
  }
  \item{fdrth}{
     \code{fdrth} is the fdr cutoff to identify statistically significant GO terms. The default is 0.05.
  }
  \item{ppth}{
    \code{ppth} is the significant level to test whether the remaining genes of the ancestor term are enriched 
    with interesting genes after removing the genes in its significant offspring terms. The default is 0.05.
  }
  \item{pcth}{
    \code{pcth} is the significant level to test whether the frequency of interesting genes in the offspring 
    terms are significantly different from that in the ancestor term. The default is 0.05. 
  }
  \item{poth}{
    \code{poth} is the significant level to test whether the overlapping genes of one term is significantly
    different from the non-overlapping genes of the term. The default is 0.05.
  }
  \item{peth}{
    \code{peth} is the significant level to test whether the non-overlapping genes of one term is enriched with
    interesting genes. The default is 0.05.
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
   There are two types of result output of GO-function. The first type is that GO-function saves a table contained all
   statistically significant terms to a CSV file (e.g. "sigTerm.csv") in the current working folder. This table contains 
   seven columns: goid, name, refnum (the number of the reference genes in a GO term), interestnum (the number of the
   interesting genes in a GO term), pvalue, adjustp (the corrected p value by the fdr control), FinalResults. The
   "FinalResults" contains three types: (1) "Local" represents terms removed after treating for local redundancy; (2)
   "Global" represents terms removed after treating for global redundancy; (3) "Final" represents the remained terms with 
   evidence that their significance should not be simply due to the overlapping genes. GO-function also saves the
   structure of GO DAG for all statistic significant terms into a plot (e.g. "sigTerm.bmp") in the current 
   working folder. In this plot, "circle", "box" and "rectangle" represent "Local", "Global" and "Final" terms
   in the table, respectively. The different color shades represent the  adjusted p values of the terms.
}

\author{
    Jing Wang
}

\note{
  GO-function use the GO data and annotation data from Bioconductor, so the user does not need to update the data manually.
}


\examples{
       
       data(exampledata)
       sigTerm <- GOFunction(interestGenes, refGenes, organism = "org.Hs.eg.db", ontology= "BP", fdrmethod = "BY", 
       fdrth = 0.05, ppth = 0.05, pcth = 0.05, poth = 0.05, peth = 0.05, bmpSize = 2000, filename="sigTerm")
}

\keyword{ methods }
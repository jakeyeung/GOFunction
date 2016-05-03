localRedundancy <-
function(sigTerm, generalAnn, sigTermRelation, annRef, annInterest, ppth, pcth)

{
     annRef <- unique(annRef[,1])
     allRefnum <- length(annRef)
     annInterest <- unique(annInterest[,1])
     allInterestnum <- length(annInterest)
     
     sigTermRelationRe <- sigTermRelation
     sigLabel <- array(0,dim=c(nrow(sigTerm),1))
     sigTerm$Label <- sigLabel
     sigTerm$SeLabel <- sigLabel
     La <- sigTerm[,1]
     noRelationTerm <- setdiff(La,union(sigTermRelation[,1],sigTermRelation[,2]))
     sigTerm[sigTerm[,1] %in% noRelationTerm,7] <- 1
     La <- setdiff(La,noRelationTerm)
     
     while (length(La) > 0) {
         
         sigTermRelationRe <- sigTermRelationRe[sigTermRelationRe[,2] %in% La,]
         if (nrow(sigTermRelationRe) == 0)
                leafnode <- La
         else
                leafnode <- setdiff(sigTermRelationRe[,2],sigTermRelationRe[,1])

         La <- setdiff(La,leafnode)
         
             
         for (j in c(1:length(leafnode))){
              node <- leafnode[j];
              genes <- generalAnn[generalAnn[,2]==node,1]
              genes <- intersect(genes,annRef)
              sgenes <- intersect(genes,annInterest)
              childnode <- sigTermRelation[sigTermRelation[,1]==node,2]
#cat("Panode:",node,"\n")
#cat("Chnode:",childnode,"\n")
              if (length(childnode)==0) {
                  sigTerm[sigTerm[,1]==node,7] <- 1
              }
              else {         
                  activeChild <- sigTerm[(sigTerm[,1] %in% childnode) & sigTerm[,7]==1,1]
                  if (length(activeChild)==0) {
                      sigTerm[sigTerm[,1]==node,7] <- 1
                  }
                  else {
                      allcgenes <- generalAnn[generalAnn[,2]==activeChild[1],1]
                      if (length(activeChild)>1) {
                          for (k in c(2:length(activeChild))) {
                              cgenes <- generalAnn[generalAnn[,2]==activeChild[k],1]
                              allcgenes <- union(allcgenes,cgenes)
                          }
                      }
                      allcgenes <- intersect(allcgenes,annRef)
                      allcsiggenes <- intersect(allcgenes,annInterest)
					  
					  if(length(setdiff(allcgenes,genes))>0){
# cat("Parent Node:",node," annotated with ",length(genes)," genes while child Nodes:",as.character(activeChild), "annotated with ",length(allcgenes), " genes.\n")
#cat("The number of different genes between parent node and child nodes is", length(setdiff(allcgenes,genes)),".\n")
#cat("This may be the problem from annoation package GOFunction used!\n")
					  }else{
					  
						  extragenes <- setdiff(genes,allcgenes)
						  extrasiggenes <- intersect(extragenes,annInterest)
						  if (length(extragenes)!=0) {
							  fp <- length(extrasiggenes)/length(extragenes)
							  fc <- length(allcsiggenes)/length(allcgenes)
#cat("FP:",fp," FC:",fc," length(extrasiggenes):",length(extrasiggenes)," length(extragenes):",length(extragenes)," length(allcsiggenes):", length(allcsiggenes), " length(sgenes):", length(sgenes), " length(genes):",length(genes)," length(allcgenes):",length(allcgenes),"\n\n\n\n")
						  
							  p <- 1-phyper(length(extrasiggenes)-1,allInterestnum,allRefnum-allInterestnum,length(extragenes),lower.tail = TRUE,log.p= FALSE)
							  pc <- 1-phyper(length(allcsiggenes)-1,length(sgenes),length(genes)-length(sgenes),length(allcgenes),lower.tail = TRUE,log.p= FALSE)
							  if (fp>=fc | p<=ppth) {
								  sigTerm[sigTerm[,1]==node,7] <- 1
								  if (pc>pcth)
										sigTerm[sigTerm[,1] %in% activeChild,8] <- sigTerm[sigTerm[,1] %in% activeChild,8]+1
							  }
						  }
					  }
                  }
             }
         }
      }
      sigTerm[,7] <- sigTerm[,7]+sigTerm[,8]
      sigTerm[sigTerm[,7]>1,7] <- 0;
      sigTermRedun <- sigTerm[sigTerm[,7]==1,c(1:6)]
      return(sigTermRedun)
}


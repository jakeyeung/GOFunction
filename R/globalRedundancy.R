globalRedundancy <-
function(generalAnn, sigTermRelation, annRef, annInterest, sigTermRedun, poth, peth)
{
     annRef <- unique(annRef[,1])
     allRefnum <- length(annRef)
     annInterest <- unique(annInterest[,1])
     allInterestnum <- length(annInterest)
     sigTermRedun$overlap = array(0,dim=c(nrow(sigTermRedun),1));
	 sigTermenv <- new.env(hash=T,parent=emptyenv())
	 assign("sigTerm",sigTermRedun,envir=sigTermenv)
     
     calculateEachTerm <- function (term1) {
		 sigTermRedun <- get("sigTerm",sigTermenv)
          gene1 <- generalAnn[generalAnn[,2]==term1,1]
          gene1 <- intersect(gene1, annRef)
          siggene1 <- intersect(gene1, annInterest)
          extrterm <- setdiff(sigTermRedun[,1], term1);
          calculateExtraTerm <- function(term2) {
               gene2 <- generalAnn[generalAnn[,2]==term2,1];
               gene2 <- intersect(gene2, annRef);
               siggene2 <- intersect(gene2, annInterest)
               po <- sigTermRelation[(sigTermRelation[,1]==term1 & sigTermRelation[,2]==term2) | (sigTermRelation[,1]==term2 & sigTermRelation[,2]==term1),]
               if (nrow(po)==0){
                    refov <- intersect(gene1,gene2);
                    if (length(refov)>0) {
                         sigov <- intersect(siggene1,siggene2)
                         extra1 <- setdiff(gene1,refov)
                         extrasig1 <- intersect(extra1, annInterest)
                         extra2 <- setdiff(gene2,refov)
                         extrasig2 <- intersect(extra2, annInterest)
                         if(length(extra2)==0){
                               return(0)
                         }
                         else{
                               pex2 <- 1-phyper(length(extrasig2)-1,allInterestnum,allRefnum-allInterestnum,length(extra2),lower.tail = TRUE,log.p= FALSE)
                               po2 <- 1-phyper(length(sigov)-1,length(siggene2),length(gene2)-length(siggene2),length(refov),lower.tail = TRUE,log.p= FALSE)
                               if(length(extra1)==0){
                                    if ((po2>poth) | (po2<=poth & pex2<=peth)){
                                        sigTermRedun[sigTermRedun[,1]==term1,7] <- 1
                                        assign("sigTerm",sigTermRedun,envir=sigTermenv)
                                     }
                                }
                                else{            
                                     pex1 <- 1-phyper(length(extrasig1)-1,allInterestnum,allRefnum-allInterestnum,length(extra1),lower.tail = TRUE,log.p= FALSE)
                                     po1 <- 1-phyper(length(sigov)-1,length(siggene1),length(gene1)-length(siggene1),length(refov),lower.tail = TRUE,log.p= FALSE)
                                     if((po1<=poth) & (pex1>peth)){
                                         if ((po2>poth) | (po2<=poth & pex2<=peth)){
                                              sigTermRedun[sigTermRedun[,1]==term1,7] <- 1
                                               assign("sigTerm",sigTermRedun,envir=sigTermenv)
                                         }
                                      }
                                 }
                        }
                   }
			   }
           }
		  lapply(extrterm,calculateExtraTerm)
     }
	 lapply(sigTermRedun[,1],calculateEachTerm)
	 sigTermRedun <- get("sigTerm",sigTermenv)
     sigTermRedun <- sigTermRedun[sigTermRedun[,7]==0,c(1:6)]
     return(sigTermRedun);
}


enrichmentFunction <-
function(annRef, annInterest, method, fdrth)

{
    allRefnum <- length(unique(annRef[,1]))
    allInterestnum <- length(unique(annInterest[,1]))
    
    allAnnterm <- unique(annRef[,2])
    allAnntermL <- length(allAnnterm)


    refCount <- tapply(annRef[,1],annRef[,2],length)
    refTerm <- levels(factor(annRef[,2]))
    refTermCount <- data.frame(goid=refTerm,refnum=array(refCount,length(refCount)))
    interestCount <- tapply(annInterest[,1],annInterest[,2],length)
    interestTerm <- levels(factor(annInterest[,2]))
    interestTermCount <- data.frame(goid=interestTerm,interestnum=array(interestCount,length(interestCount)))

    ref_interest_TermCount <- refTermCount;
    ref_interest_TermCount$interestnum = array(0,dim=c(length(ref_interest_TermCount$goid),1))
    ref_interest_TermCount[ref_interest_TermCount$goid %in% interestTermCount[,1],3]=interestTermCount$interestnum;


    n <- nrow(ref_interest_TermCount)
    pv <- array(0,dim=c(n,1))
    for (i in c(1:n)){
        p <- 1-phyper(ref_interest_TermCount[i,3]-1,allInterestnum,allRefnum-allInterestnum,ref_interest_TermCount[i,2],lower.tail = TRUE,log.p= FALSE)
        pv[i,1] <- p
    }
    ref_interest_TermCount$pvalue <- pv

    ref_interest_TermCount <- ref_interest_TermCount[order(ref_interest_TermCount[,4]),]

    if (length(grep(method,"bonferroni"))){
        adjustp <- p.adjust(ref_interest_TermCount$pvalue,method ="bonferroni");
        ref_interest_TermCount$adjustp <- adjustp
        sigTerm <- ref_interest_TermCount[ref_interest_TermCount[,5]<=fdrth,];
    }
    else {
        if (length(grep(method,"BH"))){
                adjustp <- p.adjust(ref_interest_TermCount$pvalue,method ="BH");
                ref_interest_TermCount$adjustp <- adjustp
                sigTerm <- ref_interest_TermCount[ref_interest_TermCount[,5]<=fdrth,];
        }
        else {
             if (length(grep(method,"BY"))){
                    adjustp <- p.adjust(ref_interest_TermCount$pvalue,method ="BY");
                    ref_interest_TermCount$adjustp <- adjustp
                    sigTerm <- ref_interest_TermCount[ref_interest_TermCount[,5]<=fdrth,];
             }
        }
    }
    rownames(sigTerm) <- NULL
    return(list(sigTerm=sigTerm,allTerm=ref_interest_TermCount))       
}


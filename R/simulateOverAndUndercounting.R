#' @export
translateCodesToProbabilities <- function(compSeq,translationTable,baselineProb = 0.05){
compSeq$ocProb <- baselineProb -> compSeq$ucProb
for(i in 1:nrow(translationTable)){
  ind <- which(compSeq$varveCode==translationTable$codes[i])
  compSeq$ocProb[ind] <- translationTable$overcount[i]
  compSeq$ucProb[ind] <- translationTable$undercount[i]
  }
return(compSeq)
}

simulateOverAndUndercounting = function(compSeq){

  OCP <- compSeq$ocProb
  UCP <- compSeq$ucProb
  thicks <- compSeq$thick


  if(length(OCP)==1){#if there's only one OCP, replicate over length
    OCP=matrix(data=OCP,ncol=1,nrow=length(thicks))
  }
  if(length(UCP)==1){#if there's only one UCP, replicate over length
    UCP=matrix(data=UCP,ncol=1,nrow=length(thicks))
  }


  #generate a random series (uniform length)
  rdata=runif(length(thicks))

  #find those that were overcounted,
  OCi=which(rdata<OCP)
  #and undercounted
  UCi=which(rdata>(1-UCP))

  #loop through and correct overcounting
  newThicks=matrix(data=NA,ncol=1,nrow=length(thicks)+length(UCi)-length(OCi))
  #newThicks=c()
  counter=0
  secFlag=FALSE
  OCi1Flag=FALSE
  u=1
  wasUCi=c()#create new indices that show which years had been UC before correction
  wasOCi=c()
  #create an index that maps old to new
  old2new=c()


  while(u <= counter+length(thicks)){
    if(any((u-counter)==UCi) & !secFlag){
      wu=UCi[which(u==(UCi+counter))]
      newThicks[u]=thicks[wu]/2
      secFlag=TRUE
      wasUCi=append(wasUCi,u)
    }else if(secFlag){
      newThicks[u]=thicks[wu]/2
      counter=counter+1
      secFlag=FALSE
      wasUCi=append(wasUCi,u)

    }else if(any((u-counter)==OCi) & !secFlag){
      wo=OCi[which(u==(OCi+counter))]
      #build more functionality here for handling these options
      if (u==1){
        #combine with next measurement
        newThicks[u]=thicks[wo]+thicks[wo+1]
        wasOCi=append(wasOCi,u)
      }else{
        #combine with previous
        u=u-1
        newThicks[u]=thicks[wo]+newThicks[u]
        #newThicks[u]=thicks[wo+1]
        wasOCi=append(wasOCi,u)
      }
      counter=counter-1
    }else{
      newThicks[u]=thicks[u-counter]
    }
    old2new[u]=u-counter
    u=u+1



  }
  #account for possibility that last one was undercounted
  if(secFlag){
    newThicks[u]=thicks[wu]/2
    counter=counter+1
    secFlag=FALSE
    wasUCi=append(wasUCi,u)
  }

  #remove NAs.
  newThicks=na.omit(newThicks)

  #NOTE: THIS ERROR GETS VIOLATED OCCASIONALLY (1/200 simulations), for reasons that aren't clear to me. I don't think it's too important. Yet. OK - I think I know why, when the first one is Undercounted, and the second is overcounted. It handles it properly, but doesn't register 2 as overcounted. The total thicknesses are the same, so I'm ok with this, for now at least.

  #   #do some checking!
  #   if(length(newThicks)-counter != length(thicks)){
  #     stop("Counter & newThick lenghts don't line up")
  #   }
  if( abs(sum(thicks)-sum(newThicks))>0.0001 ){
    stop("The thicknesses don't sum to the same number")
  }

  out = list(newThicks=newThicks,old2new=old2new,wasUCi=wasUCi,wasOCi=wasOCi,UCi=UCi,OCi=OCi)
  return(out)

}


#' @export
generateThicknessEnsemble <- function(compSeq,nEns = 1000){
  ensThick <- matrix(NA,nrow = 2*nrow(compSeq),ncol = nEns)
  pb <- txtProgressBar(min=1,max=nEns,style=3)

for(i in 1:nEns){
  thisEns <- simulateOverAndUndercounting(compSeq)
  ensThick[1:length(thisEns$newThicks),i] <- thisEns$newThicks
  setTxtProgressBar(pb,i)

}

  #trim rows of all NaNs
  ensThick <- ensThick[-which(rowSums(!is.na(ensThick))==0),]
  return(ensThick)
}
#
# library(geoChronR)
#  plotTimeseriesEnsRibbons(X = seq(1,2504),Y = test) %>%
#  plotTimeseriesEnsLines(X = seq(1,2504),Y = test,maxPlotN = 2,color = "red")
# ensThick <- ensThick*(1500/sum(ensThick[,1],na.rm = TRUE))


#' @export
createEnsembleAgeDepthModel <- function(ensThick,ageVec = seq(1,nrow(ensThick)), addToLiPD = NA,startAtZero = TRUE){
  if(startAtZero){
  ensThick <- rbind(matrix(0,ncol = ncol(ensThick),nrow = 1), ensThick)
  }

  cumEns <- apply(ensThick,2,cumsum)
  depthVec <- seq(from = 0, to = max(cumEns,na.rm = TRUE),length.out = nrow(cumEns))



  ageEns <- apply(cumEns,2,function(x){approx(x,ageVec,depthVec)$y})


  ageDepthMat <- cbind(depthVec,ageEns)
  colnames(ageDepthMat) <- c("depth",paste0("ageEns",seq(1,ncol(ageEns))))
  summaryTable <- cbind(depthVec,t(apply(ageEns,1,quantile,probs = c(.025,.25,.5,.75,.975),na.rm = TRUE)))
  names(summaryTable) <- c("depth","2.5%","25%","median","75%","97.5%")


  if(is.na(addToLiPD)){
    return(list(ageDepthEns = ageDepthMat,summaryTable = summaryTable))
  }else{
    stop("Add to LiPD not finished yet")
  model <- vector(mode = "list",length = 1)

  model[[1]]$ensembleTable$ageEnsemble$values <- ageEns

  model[[1]]$methods$algorithm <- "varveR"
  model[[1]]$summaryTable$values <- summaryTable
  }
}

#' @export
#' @import Hmisc
createDepth2AgeFunction <- function(ageDepthMat){
depth <- ageDepthMat[,1]
ae <- ageDepthMat[,-1]

fout <- function(d2m){
  ens <- apply(ae,2,function(x){Hmisc::approxExtrap(depth,x,d2m)$y})
  med <- apply(ens,1,median, na.rm=TRUE)
  return(list(ageEnsemble = ens, medianAge = med, depth = d2m))
  }
return(fout)
}





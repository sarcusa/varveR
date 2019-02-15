#' @export
#' @author Nick McKay
#' @title Read varve thicknesses from shape file
#' @description Pulls varve thicknesses and error codes from a shapefile counted in GIS
#' @importFrom sf st_read
#' @import dplyr tibble magrittr
#' @param filename the path to a shape file
#' @return A data.frame with the varve data
readVarveShapefile <- function(filename,codeCol = c("conf","VarveID"),markerCol = c("Marker","Markers"),varveTop = "top",scaleToThickness = NA){
  shape <-  sf::st_read(filename,quiet = TRUE)

  varveCoords <- tibble::tibble(x1 = unlist(lapply(shape$geometry,"[[",1)),
                            x2 = unlist(lapply(shape$geometry,"[[",2)),
                            y1 = unlist(lapply(shape$geometry,"[[",3)),
                            y2 = unlist(lapply(shape$geometry,"[[",4)))


  #pull out varve code metadata
  hasCode <- FALSE
  if(any(names(shape) %in% codeCol)){
    gc <- which(names(shape) %in% codeCol)
    varveCoords$varveCode = shape[[gc]]
    hasCode <- TRUE
  }else{
    warning(paste(filename,"- no varve codes. Other names include:"))
    print(names(shape))
  }

  #pull out marker layers
  hasMarker <- FALSE
  if(any(names(shape) %in% markerCol)){
    gc <- which(names(shape) %in% markerCol)
    varveCoords$markers = shape[[gc[1]]]
    hasMarker <- TRUE
  }else{
    warning(paste(filename,"- no marker layers. Other names include:"))
    print(names(shape))
  }

  if(varveTop == "top"){
    varves <- varveCoords %>%
      dplyr::arrange(desc(y1)) %>%
      dplyr::mutate(thick = abs(y1-y2)) %>%
      dplyr::filter(thick > 0) %>% #remove layers of zero thickness
      dplyr::mutate(count = seq_along(y1))
  }else if(varveTop == "right"){
    varves <- varveCoords %>%
      dplyr::arrange(desc(x1)) %>%
      dplyr::mutate(thick = abs(x1-x2)) %>%
      dplyr::filter(thick > 0) %>% #remove layers of zero thickness
      dplyr::mutate(count = seq_along(x1))
  }else if(varveTop == "left"){
    varves <- varveCoords %>%
      dplyr::arrange(x1) %>%
      dplyr::mutate(thick = abs(x1-x2)) %>%
      dplyr::filter(thick > 0) %>% #remove layers of zero thickness
      dplyr::mutate(count = seq_along(x1))
  }else if(varveTop == "bottom"){
    varves <- varveCoords %>%
      dplyr::arrange(y1) %>%
      dplyr::mutate(thick = abs(y1-y2)) %>%
      dplyr::filter(thick > 0) %>% #remove layers of zero thickness
      dplyr::mutate(count = seq_along(y1))
}

  if(hasCode & hasMarker){
    out <- dplyr::select(varves,count,thick,varveCode, markers)
  }else if(!hasCode & hasMarker){
    out <- dplyr::select(varves,count,thick, markers)
  }else if(hasCode & !hasMarker){
    out <- dplyr::select(varves,count,thick,varveCode)
  }else{
    out <- dplyr::select(varves,count,thick)
  }


  #Scale measurements to match section thickness?
  if(!is.na(scaleToThickness)){
    out$thick <- out$thick*(scaleToThickness/sum(out$thick))
  }

  return(out)
}

readVarveDir <- function(directory = NULL,codeCol = c("conf","VarveID"),varveTop = "top",markerCol = c("Marker","Markers"),scaleToThickness = NA){
if(is.null(directory)){
  print("Select a file in the desired directory")
  directory <- dirname(file.choose())
}

  f2g <- list.files(directory,pattern = "*.shp")
  vOut <- vector(mode = "list",length = length(f2g))
  for(i in 1:length(f2g)){
    vOut[[i]] <- readVarveShapefile(file.path(directory,f2g[i]),codeCol = codeCol,varveTop = varveTop,markerCol = markerCol,scaleToThickness = scaleToThickness)

  }
  names(vOut) <- f2g
return(vOut)
}


determineMarkerLayerOrder <- function(sequence){

  #take a guess at top and bottom of sequence:

  #pull out marker layers for each
  mls <- sapply(sapply(sequence,"[[","markers"),levels)

  #which ones might be top and bottom?
  tb <- which(sapply(mls,length)==1)
  top <- c()
  bot <- c()
  for(i in tb){
    #find a match
    mms <- unlist(sapply(mls,function(x){which(x==mls[[i]])}))
    if(max(mms)==1){
    top <- c(top,i)
    }
    if(max(mms)>1){
      bot <- c(bot,i)
    }
  }

  #do we have top and bottom?
  if(length(top)==1 & length(bot)==1){#YAY!
    secOrder <- rep(NA,length(sequence))
    secOrder[1] <- top
    secOrder[length(sequence)] <- bot
    secRemaining <- dplyr::setdiff(seq_along(sequence),c(top,bot))
  }else{
    stop("cant identify a top and bottom section")
  }

  #loop through and figure out the rest of the order...

  while(length(secRemaining)>0){
    remSpots <- which(is.na(secOrder))
    i <- min(remSpots)
    mms <- sapply(mls,function(x){ifelse(length(which(mls[[secOrder[i-1]]] %in% x))>0,which(mls[[secOrder[i-1]]] %in% x),0)})
    hasMatch <- setdiff(which(mms>0),secOrder)
    if(length(hasMatch)==1){#only one new match!
      secOrder[i] <- hasMatch
    }else if(length(hasMatch)>1){
      stop("multiple matches of marker layers. Should be a pair.")
    }else{
      stop("cannot find a matching marker layer")
    }
    secRemaining <- dplyr::setdiff(seq_along(sequence),na.omit(secOrder))
  }

return(secOrder)

}

combineSectionsByMarkerLayer <- function(sortedSequence){
#presently tries to handle multiple marker layer options stochastically.


#build a giant, concatenated dataframe
  allCounts <- dplyr::bind_rows(sortedSequence, .id = "section")
  sections <- unique(allCounts$section)

#get all the markerlayer names.
  mln <- na.omit(unique(allCounts$markers))

  #see if there are multiple marker layers that connect sections
  nMark <- allCounts %>%
    group_by(section) %>%
    summarize(markerLayers = sum(!is.na(markers))) %>%
    arrange(match(section,sections))

  nMark$expected <- 2
  nMark$expected[c(1,nrow(nMark))] <- 1
  #find sections with "extra" markers
  extra <- which(nMark$markerLayers-nMark$expected > 0)

  if(length(extra)>0){
    #look for similar pairs.
    extraCounts <- dplyr::bind_rows(sortedSequence[extra], .id = "section") %>%
      group_by(markers) %>%
      summarize(nMark = n())

    #summarize into a 2 column matrix, only one is needed.
    extraMarks <- matrix(na.omit(extraCounts$markers[extraCounts$nMark>1]),ncol = 2,byrow = T)


    for(i in 1:nrow(extraMarks)){
      mln <- setdiff(mln,sample(extraMarks[i,],1)) #remove one at random.
    }
  }

#start at the top to the first markerlayer
  sp <- min(which(allCounts$markers==mln[1]))
  combSeq <- allCounts[1:sp,]

  #Deal with intermediate sections
  for(n in 2:length(mln)){
    #find second instance of previous markerlayer,
    st <- max(which(allCounts$markers==mln[n-1]))

    #find first instance of next marker layer.
    sp <- min(which(allCounts$markers==mln[n]))

    if(sp <= st){#that's bad, markers out of order
      stop(paste("marker layers",mln[n-1],"&",mln[n],"out of order"))
      }

    combSeq <- bind_rows(combSeq,allCounts[(st+1):sp,])
  }


  #Then append the bottom section
  st <- max(which(allCounts$markers==mln[n]))
  combSeq <- bind_rows(combSeq,allCounts[(st+1):nrow(allCounts),])

#rename and extend tibble
combSeq <- combSeq %>%
  mutate(sectionCounts = count) %>%
  mutate(count = seq_along(count))


plotCompositeSequence(combSeq)
return(combSeq)



}

plotCompositeSequence <- function(compSeq){
  colourCount = length(unique(compSeq$section))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

  return(ggplot(compSeq)+geom_line(aes(x = 1950-count, y = thick,color=section))+theme_bw()+scale_colour_manual(values = getPalette(colourCount)))
}



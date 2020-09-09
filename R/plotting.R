#' Plot a bunch of varve sections
#'
#' @param varveSequenceList
#' @param output A pdf file to concatenate the output into
#' @import purrr pdftools
#'
#' @return a
#' @export
plotSections <- function(varveSequenceList,output = NA){
 
 for(i in 1:length(varveSequenceList)){
    
    if(anyNA(varveSequenceList[[i]]$varveCode)){
      print(paste0("check slide ", names(varveSequenceList[i])," because varveID looks wrong around varve ", which(is.na(varveSequenceList[[i]]$varveCode))))
    }
    
  }
 
 allPlots <- purrr::map2(varveSequenceList,names(varveSequenceList),plotVarveSection)
 if(!is.na(output)){
   purrr::map2(paste0(file.path(tempdir(),names(varveSequenceList)),".pdf"),allPlots,ggsave)
   apdf <- list.files(tempdir(),pattern = "*.shp.pdf",full.names = T)
   pdftools::pdf_combine(apdf,output = output)
 }
 return(allPlots)
}


#' Plot varve section
#'
#' @param section
#' @param sectionName
#' @param varveCodeScale
#' @import ggplot2 dplyr pracma
#' @return a ggplot list
#' @export
#'
#' @examples
plotVarveSection <- function(section,sectionName, varveCodeScale = 0.8){
  #scale varve codes to counts
  section$scaledVarveCode <-  varveCodeScale*section$varveCode*max(section$thick)/max(section$varveCode)
  svcLabels <- c(as.character(unique(sort(section$varveCode))),"Varve Codes")
  svcLabelY <- varveCodeScale*unique(sort(section$varveCode))*max(section$thick)/max(section$varveCode)
  svcLabelY <- c(svcLabelY,max(svcLabelY)*1.1)
  svcLabelX <- quantile(section$count,probs = 1)

# nearest neighbor varve codes for plotting.
  iCount <- seq(from = min(section$count),to = max(section$count),by = 0.1)
  nn <- data.frame(count = iCount,
                   varveCode = pracma::interp1(section$count,section$scaledVarveCode,xi = iCount,method = "nearest" ))


  #get marker layers for plotting
  ml <- filter(section, !is.na(markers))


  secPlot <- ggplot()+
    #plot varve codes
    geom_area(data = nn, aes(x = count, y = varveCode),fill = "Gray",alpha = 0.5)+
    #add varve code scale
    geom_label(aes(x= svcLabelX,y = svcLabelY, label = svcLabels),colour = "Gray")+
    geom_line(data = section, aes(x = count, y = thick))+
    ggtitle(sectionName)+
    ylab("Varve Thickness (arbitary units)")+
    xlab("Count")+
    theme_bw()

  if(nrow(ml)){
    secPlot <- secPlot+
    #add marker layers
    geom_vline(data = ml, aes(xintercept = count),color = "red")+
    geom_label(data = ml, aes(x= count,y = 0, label = markers),colour = "Red")
  }


  return(secPlot)


}

#' plot a composite sequence
#'
#' @param compSeq
#' @import RColorBrewer
#' @return a plot
#' @export
plotCompositeSequence <- function(compSeq){
  colourCount = length(unique(compSeq$section))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

  return(ggplot(compSeq)+geom_line(aes(x = count, y = thick,color=section))+theme_bw()+scale_colour_manual(values = getPalette(colourCount)))
}


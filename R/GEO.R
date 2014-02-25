###http://www.ncbi.nlm.nih.gov/geo/info/geo_illu.html
createGEOMeta <- function(summaryData){
  
  data(GEOmetaTemplate)
  
  pd <- pData(summaryData)
  metanames <- metaTemplate[["SAMPLES"]]
  
  metaMatrix <- matrix(nrow = nrow(pd),ncol=length(metanames))
  colnames(metaMatrix) <- metanames
  for(i in 1:length(metanames)){
    if(metanames[i] %in% colnames(pd)) metaMatrix[,metanames[i]] <- pd[,metanames[i]]
  }
  
  if(all(is.na(metaMatrix[,"Sample name"]))) metaMatrix[,"Sample name"] <- sampleNames(summaryData)
 
  extraMeta <- pd[,setdiff(colnames(pd), metanames)]
  colnames(extraMeta) <- paste0("characteristics:",colnames(extraMeta))
  
  metaMatrix <- cbind(metaMatrix, extraMeta)
  
  newMeta <- metaTemplate
  newMeta[["SAMPLES"]] <- metaMatrix

  
  
  cat("SERIES\n",file="GEOmeta.txt")
  sapply(newMeta["SERIES"], function(x) cat(paste0(as.character(x), "\n"), file="GEOmeta.txt",append=T))
  cat("SAMPLES\n",file="GEOmeta.txt",append=T)
  cat(paste0(paste(colnames(metaMatrix),collapse="\t"),"\n"),file="GEOmeta.txt",append=T)
  apply(metaMatrix, 1, function(x) cat(paste(x, collapse="\t"),"\n",file="GEOmeta.txt",append=T))
  cat("PROTOCOLS\n",file="GEOmeta.txt",append=T)
  sapply(newMeta["PROTOCOLS"], function(x) cat(paste0(as.character(x), "\n"), file="GEOmeta.txt",append=T))
  
  
  
}

createGEOMatrix <- function(data,forceDetection=TRUE){
  
  annoPkg <- paste("illumina", annotation(data), ".db",sep="")
  
  
  annoLoaded <- require(annoPkg, character.only=TRUE)
  
    if(annoLoaded){
    
      mapEnv <-  as.name(paste("illumina", annotation(summaryData), "ARRAYADDRESS",sep=""))
    
      features <- mappedkeys(eval(mapEnv))
    
  
      summaryData <- data[features,]
    
      if(dim(summaryData)[1] == 0) message("Annotation map not successful")
      
      eMat <- exprs(summaryData)
      
      detMat <- matrix(nrow=nrow(eMat),ncol=ncol(eMat),NA)
      
      det <- Detection(summaryData)
      
      if(is.null(det)){
            
        if(forceDetection) {
          message("Calculating Detection scores")
          det <- try(calculateDetection(summaryData))
    
          if(class(det ) == "try-error") {
            message("Could not calculate detection scores")
          }
          else detMat <- det
        }
            
        } else detMat <- det
      
      
      fullMat <- matrix(nrow=nrow(eMat), ncol=ncol(eMat)*2)
      fullMat[,seq(1, by = 2, len = ncol(eMat))] <- eMat
      fullMat[,seq(2, by = 2, len = ncol(eMat))] <- det
      colnames(fullMat)[seq(2, by = 2, len = ncol(eMat))] <- "Detection Pval"
      
      colnames(fullMat)[seq(1, by = 2, len = ncol(eMat))] <- colnames(eMat)
      fullMat <- cbind("ID_REF"= features, fullMat)
      
      write.table(fullMat, file="GEOProcessedMatrix.txt",sep="\t",quote=F,row.names=F)
      
    
  } else message("Could not load required annotation package ",annoPkg)
  
}
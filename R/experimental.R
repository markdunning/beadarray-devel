

simpleDE <- function(summaryData,SampleGroup=NULL,...){
  
  if(is.null(SampleGroup)) SampleGroup <- pData(summaryData)$SampleGroup
  
  design <- model.matrix(~0+as.factor(SampleGroup))
  colnames(design) <- as.character(levels(as.factor(SampleGroup)))
  
  contrast <- vector()
  
  for (a in 1:length(levels(as.factor(SampleGroup)))){
    for (b in 1:length(levels(as.factor(SampleGroup)))){ 
      if (a!=b){
        if (a<b){
          contrast[length(contrast)+1] <- paste(levels(as.factor(SampleGroup))[a],levels(as.factor(SampleGroup))[b],sep="-")
        }
      }
    }
  }
  
  fit <- lmFit(exprs(summaryData), design)
  cnt <- paste(colnames(design)[1],colnames(design)[2],sep="-")
  cMat <- makeContrasts(contrasts =contrast,levels=design)
  fit2 <- contrasts.fit(fit,cMat)
  efit <- eBayes(fit2)
  
  
  summaryData@deResults <- assayDataNew(LogFC = efit$coef, PValue = efit$p.value, LogOdds = efit$lods)
  summaryData
  
}



setGeneric("SampleGroup", function(object) standardGeneric("SampleGroup"))

setMethod("SampleGroup", signature(object = "ExpressionSetIllumina"), function(object) object@SampleGroup)
setGeneric("SampleGroup<-", function(object, value) standardGeneric("SampleGroup<-"))

setReplaceMethod("SampleGroup",
                 signature=signature(
                   object="ExpressionSetIllumina",
                   value="character"),
                 function(object, value) {
                   object@SampleGroup <- value
                   object
                 })

ExpressionSetFromGEO <- function(gse){
  
summaryData <- new("ExpressionSetIllumina")
exprs(summaryData) <- exprs(gse)
phenoData(summaryData) <- phenoData(gse)
summaryData@channelData[[1]] <- rep("G", length(sampleNames(gse)))
featureData(summaryData) <- featureData(gse)[,1:3]

annotation(summaryData) <- switch(annotation(gse), GPL6947="Humanv3", GPL10558="Humanv4", GPL6887="Mousev2", GPL6102="Humanv2")
summaryData <- addFeatureData(summaryData)

summaryData
}

setSampleGroup <- function(summaryData, SampleGroup){
  pData(summaryData)$SampleGroup <- SampleGroup
  summaryData
}

collapse <- function(summaryData, fn = function(x) IQR(x,na.rm=T), unit = "SYMBOL",stats=NULL,units=NULL){
  
  fd <- fData(summaryData)
  unitCol <- grep(unit, colnames(fd))
  units <- na.omit(unique(fd[,unitCol]))
  
  if(is.null(stats)) stats <- apply(exprs(summaryData), 1, fn)
  
  
  probeOrder <- fd[order(stats,decreasing=T),]
  
  chosenID <- rownames(probeOrder[match(units, probeOrder[,unitCol]),])
  
  collapsedData <- summaryData[match(chosenID,featureNames(summaryData)),]
}


toRangedData <- function(summaryData){
  
  annoName <- annotation(summaryData)
  
  annoLoaded <- require(paste("illumina", annoName, ".db",sep=""), character.only=TRUE)
  
   if(annoLoaded){
    
    
      mapEnv <-  as.name(paste("illumina", annoName, "GENOMICLOCATION",sep=""))
      fn <- featureNames(summaryData)
      fn <- fn[which(fn %in% mappedkeys(eval(mapEnv)))]
        
      locs <- mget(fn,eval(mapEnv),ifnotfound=NA)
      
      locs <- lapply(locs, function(x) gsub(" ", ",", x,fixed=T))
      
      asLocMatrix <- function(str){
        x<- do.call("rbind",sapply(strsplit(as.character(str), ",",fixed=T)[[1]], function(x) as.vector(strsplit(x, ":",fixed=T))))
      }
    
      locMat <- lapply(locs, asLocMatrix)
      
      rn <- rep(names(locs), unlist(lapply(locMat, nrow)))
      
      locMat <- do.call("rbind", locMat)
      rng <- GRanges(locMat[,1], IRanges(as.numeric(locMat[,2]), as.numeric(locMat[,3]),names=rn),strand=locMat[,4])
      values(rng) <- fData(summaryData)[match(rn, featureNames(summaryData)),]
      rng
  }
  
}

volcanoplot <- function(summaryData){
  
  data <- data.frame(LogFC = LogFC(summaryData), LogOdds = LogOdds(summaryData), fData(summaryData))
  colnames(data)[1] <- "LogFC"
  colnames(data)[2] <- "LogOdds"
  
  ggplot(data, aes(x = LogFC, y = LogOdds)) + geom_point(color="steelblue",alpha=0.3)
  
}

plotProbe <- function(symbol,rng, tx){
  data(genesymbol)  
  p1 <- autoplot(tx, which=genesymbol[symbol])
  p2 <- rng[which(rng %over% genesymbol[symbol])]
  tracks(p1, autoplot(p2,aes(fill=PROBEQUALITY)))
  
}

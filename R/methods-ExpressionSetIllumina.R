
setMethod("initialize", "ExpressionSetIllumina",
          function(.Object,
                   assayData = assayDataNew(exprs=exprs,se.exprs=se.exprs, nObservations=nObservations, Detection=Detection, storage.mode="list"),

                   phenoData = new("AnnotatedDataFrame"),
                   exprs=new("matrix"),
                   se.exprs=new("matrix"),
                   nObservations=new("matrix"),
                   Detection=new("matrix"),
                   annotation = character(),
                   featureData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME")
  )
 {
            .Object<-callNextMethod(.Object,
                           assayData = assayData,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation,
                           featureData = featureData
			   )
            
            .Object
          })




setGeneric("LogFC", function(object) standardGeneric("LogFC"))

setMethod("LogFC", signature(object="ExpressionSetIllumina"), function(object) assayDataElement(object@deResults, "LogFC"))

setGeneric("LogFC<-", function(object, value) standardGeneric("LogFC<-"))


setGeneric("LogOdds", function(object) standardGeneric("LogOdds"))

setMethod("LogOdds", signature(object="ExpressionSetIllumina"), function(object) assayDataElement(object@deResults, "LogOdds"))

setGeneric("LogOdds<-", function(object, value) standardGeneric("LogOdds<-"))



setMethod("[", "ExpressionSetIllumina", function(x, i, j, ..., drop = FALSE) {
  if (missing(drop))
    drop <- FALSE
  
  if (!missing(j)) {
    phenoData(x) <- phenoData(x)[j,, ..., drop = drop]
    protocolData(x) <- protocolData(x)[j,, ..., drop = drop]
  }
  
  if (!missing(i))
    featureData(x) <- featureData(x)[i,,..., drop=drop]
  ## assayData; implemented here to avoid function call
  orig <- assayData(x)
  ###I took this code from the eSet methods in Biobase to allow for empty se.exprs, nObservations, Detection
  storage.mode <- Biobase:::assayDataStorageMode(orig)
  
  assayData(x) <-
    switch(storage.mode,
           environment =,
           lockedEnvironment = {
             aData <- new.env(parent=emptyenv())
             if (missing(i))                     # j must be present
               for(nm in ls(orig)) {
                 if(nrow(orig[[nm]])>0)  aData[[nm]] <- orig[[nm]][, j, ..., drop = drop]
                 else aData[[nm]] <- orig[[nm]]
               }
             
             else {                              # j may or may not be present
               if (missing(j))
                 for(nm in ls(orig)){
                   if(nrow(orig[[nm]])>0) aData[[nm]] <- orig[[nm]][i,, ..., drop = drop]
                    else aData[[nm]] <- orig[[nm]]
                 }
                else
                 for(nm in ls(orig)){ 
                   if(nrow(orig[[nm]])>0) aData[[nm]] <-orig[[nm]][i, j, ..., drop = drop]
                   else aData[[nm]] <- orig[[nm]]
                 }
             }
             if ("lockedEnvironment" == storage.mode) assayDataEnvLock(aData)
             aData
           },
           
           list = {
             if (missing(i))                     # j must be present
               lapply(orig, function(obj) {
                 if(nrow(obj)>0) obj[, j, ..., drop = drop]
                 else obj
               })
                           
             else {                              # j may or may not be present
               if (missing(j))
                 lapply(orig, function(obj){ 
                   if(nrow(obj)>0) obj[i,, ..., drop = drop]
                   else obj
                 })
                 
               
               else
                 lapply(orig, function(obj){
                   if(nrow(obj)>0) obj[i, j, ..., drop = drop]
                   else obj
                 }) 
             }
           }
         )
  
  
  
  

	if(!is.null(x@QC) && !missing(j)) x@QC<-x@QC[j,]
    
  if(!is.null(x@channelData) && missing(j)) x@channelData<-x@channelData
	
  de <- try(x@deResults, silent=T)

  if(class(de)!= "try-error" & length(x@deResults)>0){	

  if(!missing(i)){
    # i must be present to subset the DE results
  orig <- x@deResults
  aData <- new.env(parent=emptyenv())

  	if (missing(j)){
    
       	for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i,, ..., drop = drop]
     	}
     	else {
       	for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i, j, ..., drop = drop]
  	}
          
  x@deResults <- aData
  }  
  
}
         x

})




setValidity("ExpressionSetIllumina", function(object) {
  assayDataValidMembers(assayData(object), c("exprs", "se.exprs", "nObservations"))
})

setMethod("dim", "ExpressionSetIllumina", function(x) {

	nFeatures = nrow(fData(x))
	nSamps = length(sampleNames(x))
	nChannels = length(channelNames(x))

    c("Features"=nFeatures, "Samples"=nSamps, "Channels"=nChannels)
 } )



setMethod("exprs", signature(object="ExpressionSetIllumina"), function(object) assayDataElement(object, "exprs"))

#setGeneric("exprs<-", function(object, value) standardGeneric("exprs<-"))

setReplaceMethod("exprs", signature(object="ExpressionSetIllumina",value="matrix"), function(object, value){
	assayDataElementReplace(object, "exprs", value)
})

setMethod("se.exprs", signature(object="ExpressionSetIllumina"), function(object) assayDataElement(object, "se.exprs"))

#setGeneric("se.exprs<-", function(object, value) standardGeneric("se.exprs<-"))

setReplaceMethod("se.exprs", signature(object="ExpressionSetIllumina",value="matrix"), function(object, value){
	assayDataElementReplace(object, "se.exprs", value)
})

setGeneric("nObservations", function(object) standardGeneric("nObservations"))

setMethod("nObservations", signature(object="ExpressionSetIllumina"), function(object) assayDataElement(object, "nObservations"))

setGeneric("nObservations<-", function(object, value) standardGeneric("nObservations<-"))

setReplaceMethod("nObservations", signature(object="ExpressionSetIllumina",value="matrix"), function(object, value){
	assayDataElementReplace(object, "nObservations", value)
})

setGeneric("Detection", function(object) standardGeneric("Detection"))

setMethod("Detection", signature(object="ExpressionSetIllumina"), function(object) assayDataElement(object, "Detection"))

setGeneric("Detection<-", function(object, value) standardGeneric("Detection<-"))

setReplaceMethod("Detection", signature(object="ExpressionSetIllumina",value="matrix"), function(object, value){
	assayDataElementReplace(object, "Detection", value)
})

setMethod("show", signature(object="ExpressionSetIllumina"), function(object) {
  callNextMethod(object)
  
  cat("QC Information\n")
  cat(" Available Slots:  ")
  cat(names(object@QC))
  nms=selectSome(colnames(object@QC@data))
  cat("\n  QC Items:", paste(nms, collapse=", "))
  nms=selectSome(sampleNames(object@QC))
  cat("\n  sampleNames:", paste(nms, collapse=", "))
  cat("\n")
})


setGeneric("qcData", function(object) standardGeneric("qcData"))

setMethod("qcData", signature(object="ExpressionSetIllumina"), function(object) object@QC@data)



#setGeneric("exprs<-", function(object, value) standardGeneric("exprs<-"))

#setReplaceMethod("exprs", "ExpressionSetIllumina", function(object, value){
#	assayDataElementReplace(object, "exprs", value)
#})

#setReplaceMethod("exprs", c("ExpressionSetIllumina", "matrix"), function(object, value) {
#  assayDataElementReplace(object, "exprs", value)
#})

#setReplaceMethod("se.exprs", c("ExpressionSetIllumina", "matrix"), function(object, value) {
#  assayDataElementReplace(object, "se.exprs", value)
#})

.mergeAssayData<-function(x, y, newdimnames) {
  # this is derived from assayData combine method
  # differences:
  # - allows different number of reporters/features
  # - will merge data from identical column names into 1 column ie rbind()) 
  # - only works on 2-dimensional assayData elements



  combineElement <- function(x, y) {
    outarr<-array(NA,dim=c(length(newdimnames[[1]]),length(newdimnames[[2]])),newdimnames)
    mode(outarr)<-mode(x)
    outarr[rownames(y),colnames(y)]<-y
    outarr[rownames(x),colnames(x)]<-x
    outarr
  }


  storage.mode <- storageMode(x)
  nmfunc <- assayDataElementNames

  if (storageMode(y) != storage.mode)
    stop(paste("assayData must have same storage, but are ",
               storage.mode, ", ", storageMode(y), sep=""))
  if (length(nmfunc(x)) != length(nmfunc(y)))
    stop("assayData have different numbers of elements:\n\t",
         paste(nmfunc(x), collapse=" "), "\n\t",
         paste(nmfunc(y), collapse=" "))
  if (!all(nmfunc(x) == nmfunc(y)))
    stop(paste("assayData have different element names:",
               paste(nmfunc(x), collapse=" "),
               paste(nmfunc(y), collapse=" "), sep="\n\t"))
               
  for (nm in nmfunc(x)) {
    x<-assayDataElementReplace(x,nm,combineElement(assayDataElement(x,nm),assayDataElement(y,nm)))
  }
  x
}

.mergePhenodata<-function(x , y, samples) {
  variables<-union(colnames(pData(x)),colnames(pData(y)))
  outarr<-array(data=NA,dim=c(length(samples),length(variables)),dimnames=list(samples,variables))
  outarr[sampleNames(y),colnames(pData(y))]<-as.matrix(pData(y))
  outarr[sampleNames(x),colnames(pData(x))]<-as.matrix(pData(x))
  pd<-data.frame(outarr)
  vardescs<-union(colnames(varMetadata(x)),colnames(varMetadata(y)))
  outarr<-array(data=NA,dim=c(length(variables),length(vardescs)),dimnames=list(variables,vardescs))
  outarr[colnames(pData(y)),colnames(varMetadata(y))]<-as.matrix(varMetadata(y))
  outarr[colnames(pData(x)),colnames(varMetadata(x))]<-as.matrix(varMetadata(x))
  vd<-data.frame(outarr)
  new("AnnotatedDataFrame", data=pd, varMetadata=vd)
}


#setMethod("combine", signature(x="ExpressionSetIllumina",y="ExpressionSetIllumina"), function(x, y, ...) {
setMethod("combine", signature(x="ExpressionSetIllumina",y="ExpressionSetIllumina"), function(x, y) {


  if (class(x) != class(y))
    stop(paste("objects must be the same class, but are ",
               class(x), ", ", class(y), sep=""))
  newdimnames<-list(union(featureNames(x),featureNames(y)),union(colnames(exprs(x)),colnames(exprs(y))))
  x <- .mergeAssayData(x, y, newdimnames)
  # a bit of a hack to only keep the union, and discard double entries

  newsamplenames = union(sampleNames(x), sampleNames(y))
 	
  phenoData(x) <- .mergePhenodata(x, y, newsamplenames)

  experimentData(x) <- combine(experimentData(x),experimentData(y))



    
  ## annotation -- constant
  if (any(annotation(x) != annotation(y))) {
    warning("objects have different annotations: ",
         annotation(x), ", ", annotation(y))
    annotation(x)<-unique(c(annotation(x),annotation(y)))
  }


  ##Preserve the channel names of the resulting object

  x@channelData[[1]] = c(x@channelData[[1]],y@channelData[[1]])

  x@QC = combine(x@QC,y@QC)
  
  x

})


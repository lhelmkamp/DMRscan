
#' Find differentially methylated regions using the scan statistic.
#'
#' @param meth input data, a  methylSigData object.
#' @param xvalues either "Index" to use the position values associated with the data, or "Index" to treat the data points as equally spaced along the chromosome.
#' @param dir is the temporary directory where SatScan files will be written.
#' @param smallestregion specifies the smallest region (in number of data points) that will be returned as a DMR. Default=4.
#' @param plotresult will output a graphical representation of each DMR to a pdf file in the specified directory.   Must specify build.
#' @return a list containing matrix of DMR results, as well as the likelihood values needed for plotting.  
#' @seealso \code{\link{methdiffSatScan}} wraps this function, and allows its use on data with multiple chromosomes.  
#' @export
#' @examples
#' SatScanrun(meth)
#' 
#' 
#' 
#' 

SatScanrun <- function(meth, xvalues=c("Index", "Position"), smallestregion, dir=tmpdir ) {
  
  ###################### do sitewise analysis using LRT calcs in methylSig ######################  
  myDiffSigboth <-methylSigCalc(meth) # parameters go to default unless specified in the "...".
  
  
  
  ####################### normalize LRT values via quantiles #############################
  
  logliks<-myDiffSigboth@results[,"logLikRatio"]
  dfs<-myDiffSigboth@results[,"df"]
  ts<-pt(logliks, dfs) # this gets t-distribution values from the loglikelihood and the df
  
  percentile <- ecdf(ts) # this makes a percentile function
  res<-percentile(ts) # this gets the percentile from the t values
  normlambdas<-qnorm(res) # this gives the normal quantiles.
  
  # prep the loglikelihood values for use in plotting.
  loglikmat<-data.frame(chr=myDiffSigboth@data.chr, pos=myDiffSigboth@data.start, loglik=logliks, normlambda=normlambdas)
  
  ####################### on these normalized lambdas, run RSatScan #############################
  
  # replace Infs with high value.
  #table(normlambdas==Inf)
  normlambdas[which(normlambdas==Inf)]<-10  # what do we want this to be?
  
  normlambdas.nona<-normlambdas[which(!is.na(normlambdas))]
  
  # prep for satscan
  ID<-paste(myDiffSigboth@data.chr, myDiffSigboth@data.start, sep="_")
  cases<-rep(1, length(normlambdas.nona))
  value<-normlambdas.nona
  
  # case file
  NormScanStatcas<-data.frame(ID, cases, value)
  
  # coordinates file
  if (xvalues=="Index"){
    X<-seq(from=1, to=length(normlambdas.nona))
  }
  if (xvalues=="Position"){
    X<-myDiffSigboth@data.start
  }
  Y<-rep(0, length(normlambdas.nona))
  NormScanStatgeo<-data.frame(ID, X, Y)
  
  
  # For good style, an analysis would begin by resetting the paremeter file: 
  invisible(ss.options(reset=TRUE))
  # ss.options()
  
  # Then, one would change parameters as desired
  ss.options(list(CaseFile="NormScanStat.cas", PrecisionCaseTimes=0)) # time = none
  ss.options(list(CoordinatesFile="NormScanStat.geo", AnalysisType=1, ModelType=5, TimeAggregationUnits=0))
  ss.options(list(CoordinatesType=0)) # cartesian instead of lat/long
  ss.options(list(MaxSpatialSizeInPopulationAtRisk=1)) #1% of the length can be in a zone.  (or, 1.5k sites) 
  
  # It might be reasonable at this point to check what the parameter file looks like:
  # head(ss.options(),3)
  
  
  # Then, we write the parameter file, the case file, and the geometry file to some writeable location in the OS, using the functions in package. These ensure that SaTScan-readable formats are used.
  write.ss.prm(dir, "NormScanStat")
  write.cas(NormScanStatcas, dir, "NormScanStat")
  write.geo(NormScanStatgeo, dir, "NormScanStat")
  
  
  # Now we're ready to run SaTScan. The default locations in the function correspond to the location of the executable on my disk, so those or not specified here.
  NormScanStat = satscan(dir, "NormScanStat", verbose=T) # default program location works now.
  # started at 11:25 finished at 11:26
  
  # remove the files 
  file.remove(paste(dir, "/NormScanStat.cas", sep=""))
  file.remove(paste(dir, "/NormScanStat.geo", sep=""))
  file.remove(paste(dir, "/NormScanStat.prm", sep=""))
  
  # The rsatscan package provides a summary method for satscan objects.
  summary(NormScanStat)
  NormScanStat$col
  NormScanStat$gis
  NormScanStat_results<-NormScanStat$col[which(NormScanStat$col$P_VALUE<0.05 & NormScanStat$col$NUMBER_LOC>=smallestregion),]
  cat("\nThere were", dim(NormScanStat_results)[1], "clusters with p < .05 and at least", smallestregion ,"sites.")
  pos.start<-vector()
  pos.stop<-vector()
  ind.start<-vector()
  ind.stop<-vector()
  if (dim(NormScanStat_results)[1]>0){
    for (res in 1:dim(NormScanStat_results)[1]){ 
      sites<-as.character(NormScanStat$gis[which(NormScanStat$gis$CLUSTER==NormScanStat_results$CLUSTER[res]),"LOC_ID"] )
      pos.start[res]<-min(as.numeric(gsub("^.*?_","",sites)))
      pos.stop[res]<-max(as.numeric(gsub("^.*?_","",sites)))
      inds<-which(ID %in% sites)
      
      ind.start[res]<-min(inds)
      ind.stop[res]<-max(inds)
    } # end getting sites that are in each region.
    length<-pos.stop-pos.start+1
    nsites<-ind.stop-ind.start+1
    chr<-gsub("_.*$","",NormScanStat_results$LOC_ID)
    mean_in<-NormScanStat_results[,"MEAN_IN"]
    mean_out<-NormScanStat_results[,"MEAN_OUT"]
    variance<-NormScanStat_results[,"VARIANCE"]
    pval<-NormScanStat_results[,"P_VALUE"]
    
    resultmat<-data.frame(chr, pos.start, pos.stop, ind.start, ind.stop, length, nsites, mean_in, mean_out, variance, pval )
  } else (resultmat<-vector())
  return(list(resultmat=resultmat, logliks=loglikmat))
}








#' Find differentially methylated regions using the scan statistic.
#'
#' @param meth input data, a  methylSigData object.
#' @param mydir a directory into which temporary files can be written.
#' @param xvalues defaults to "Position" to use the position values associated with the data.  Use "Index" to treat the data points as equally spaced along the chromosome.
#' @param splitcentromere splits data at the centromere within each chromosome.  May be necessary for large datasets.  Must also specify build.
#' @param build can be either hg18, hg19, or hg38.  This is necessary when splitcentromere=TRUE and/ or plotresult=TRUE.
#' @param smallestregion specifies the smallest region (in number of data points) that will be returned as a DMR. Default=4.
#' @param plotresult will output a graphical representation of each DMR to a pdf file in the specified directory.   Must specify build.
#' @return a list containing matrix of DMR results, as well as other values and parameters which can be passed to plotSatScanresult.  
#' @seealso \code{\link{SatScanrun}} which this function wraps, and \code{\link{plotSatScanresult}} which plots the results of this function. \code{\link{plotSatScanresult}} can be called within this function with plotresult=TRUE. 
#' @export
#' @examples
#' result<-DMRscan(meth=meth, mydir=getwd(), xvalues="Index", 
#' splitcentromere=TRUE, build="hg18", smallestregion=4, plotresult=TRUE)
#' 




DMRscan <- function(meth, mydir=NULL, xvalues="Position", splitcentromere=FALSE, build=NULL, smallestregion=4, plotresult=FALSE, ... ){
  
  library("rsatscan")
  library(methylSig)
  
  ##################### get the data in the right format
  
  if(class(meth)=="methylSigData"){
    cat("Data is a methylSigData object.\n")
  }else{
    stop("Data is not in a methylSigData object.  Please see the methylSigReadData function in the methylSig package.  ")
  }
  
  ##################### get/make the directory 
  if(is.null(mydir)){mydir<-getwd()}
  dir.create(paste(mydir, "/DMRScan_tmp", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  tmpdir<-paste(mydir, "/DMRScan_tmp", sep="")
  
  
  ##################### 1: splitting on centromere.
  #####################if we want to split by chromosome, see which chromosomes are present in the data
  if( splitcentromere==TRUE){
    
    
    ##################### if we want to split at centromere, read in centromeres for proper build
    # try to get build from methylSigData object if it's not specified. 
    if (is.null(build)){build<-substr(gsub("^.*?assembly=","",meth@options), 1, 4)} # this returns the first for char of the options string if no assembly is specified (eg "maxC" from maxCount=500); okay
    
    if (build %in% c("hg18", "hg19", "hg38")){
      library(GWASTools)
      if (build=="hg18"){
        data(centromeres.hg18)
        centromeres<-centromeres.hg18
      }
      if (build=="hg19"){
        data(centromeres.hg19)
        centromeres<-centromeres.hg19
      }
      if (build=="hg38"){
        data(centromeres.hg38)
        centromeres<-centromeres.hg38
      }
      
    }else{
      stop("Please specify build as hg18, hg19, or hg38 to split data at the centromere.")
    }
    
    
    found<-vector()
    loglik<-vector()
    
    mychrindlist<-rownames(table(as.character(meth@data.chr))) # this is not in order (1, 10, 11,...) instead of (1,2,..)
    mychrlistorder<-intersect(c(paste("chr",1:22, sep=""),"chrX", "chrY"), mychrindlist )
    
    # can't split by centromere if chromosomes are bad...
    badchrs<-setdiff( mychrindlist , c(paste("chr",1:22, sep=""),"chrX", "chrY") )
    goodchrs<-intersect(mychrlistorder, c(paste("chr",1:22, sep="")) )
    
    for (i in 1:length(goodchrs)){
      mychr<-goodchrs[i]
      chrnum<-gsub("chr", "", mychr)
      centromere.start<-centromeres[which(centromeres$chrom==chrnum),"left.base"]
      centromere.end<-centromeres[which(centromeres$chrom==chrnum),"right.base"]
      meth_chr<-meth[which(meth@data.chr==mychr),]
      
      for (half in 1:2){
        
        if (half==1){
          cat("\n\n\n Now on",mychr, "part", half ,"\n\n\n" )
          
          if (min(meth_chr@data.start)<=centromere.start){
            
            meth_chr_half<-meth_chr[which(meth_chr@data.start<=centromere.start),]
            result<-SatScanrun(meth=meth_chr_half, xvalues=xvalues, smallestregion=smallestregion, dir=tmpdir)
            foundchr<-result$resultmat
            loglikchr<-result$logliks
            found<-rbind(found, foundchr)
            loglik<-rbind(loglik, loglikchr)
            
          }}
        if (half==2){
          cat("\n\n\n Now on",mychr, "part", half ,"\n\n\n" )
          
          if (max(meth_chr@data.start)>=centromere.end){
            
            meth_chr_half<-meth_chr[which(meth_chr@data.start>=centromere.end),]
            result<-SatScanrun(meth=meth_chr_half, xvalues=xvalues, smallestregion=smallestregion, dir=tmpdir)
            foundchr<-result$resultmat
            loglikchr<-result$logliks
            found<-rbind(found, foundchr)
            loglik<-rbind(loglik, loglikchr)
          }
        }
        
      } # end half chrom  
    } # end chrom
    
    
    
    
    for (i in 1:length(badchrs)){
      mychr<-badchrs[i]
      chrnum<-gsub("chr", "", mychr)
      meth_chr<-meth[which(meth@data.chr==mychr),]
      
      cat("\n\n\n Now on",mychr,"\n\n\n" )
      
      result<-SatScanrun(meth=meth_chr, xvalues=xvalues, smallestregion=smallestregion, dir=tmpdir)
      foundchr<-result$resultmat
      loglikchr<-result$logliks
      found<-rbind(found, foundchr)
      loglik<-rbind(loglik, loglikchr)
      
    } # end chrom (not in 1-22, X, Y)
    
    
    
  } # end splitting at centromere
  
  
  
  
  ##################### 2: not splitting on centromere.
  
  if( splitcentromere==FALSE){
    found<-vector()
    loglik<-vector()
    
    mychrindlist<-rownames(table(as.character(meth@data.chr))) # this is not in order (1, 10, 11,...) instead of (1,2,..)
    mychrlistorder<-intersect(c(paste("chr",1:22, sep=""),"chrX", "chrY"), mychrindlist )
    
    if (length(mychrlistorder)==length(mychrindlist)){
      mychrlistfinal<-mychrlistorder
    } else {
      mychrlistfinal<-mychrindlist
    }
    
    for (i in 1:length(mychrlistfinal)){
      mychr<-mychrlistfinal[i]
      cat("\n\n\n Now on",mychr, "\n\n\n" )
      meth_chr<-meth[which(meth@data.chr==mychr),]
      
      
      result<-SatScanrun(meth=meth_chr, xvalues=xvalues, smallestregion=smallestregion, dir=tmpdir)
      foundchr<-result$resultmat
      loglikchr<-result$logliks
      found<-rbind(found, foundchr)
      loglik<-rbind(loglik, loglikchr)
      
      found<-rbind(found, foundchr)
      
    } # end chr
    
    
  } # end chr, no centromere
  
  
  
  # remove the directory too
  unlink(tmpdir,recursive = TRUE , force = TRUE)
  
  result<-list()
  result$DMRs<-found
  result$likelihoods<-loglik
  result$build<-build
  result$dir<-mydir
  
  result<<-result  # this makes result a global variable in case the plotting function fails.
  
  if (plotresult==TRUE){
    try(plotSatScanresult(result))
  }
  
  return(result)
}





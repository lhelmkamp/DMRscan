


#' Find differntially methylated regions using the scan statistic.
#'
#' @param meth input data, a  methylSigData object.
#' @param mydir a directory into which temporary files can be written.
#' @param xvalues either "Index" to use the position values associated with the data, or "Index" to treat the data points as equally spaced along the chromosome.
#' @param splitchr when true, splits the data by chromosome and analyzes each separately.
#' @param splitcentromere splits data at the centromere within each chromosome.  May be necessary for large datasets.  Must also specify build.
#' @param build can be either hg18, hg19, or hg38.  This is necessary when splitcentromere=TRUE and/ or plotresult=TRUE.
#' @param smallestregion specifies the smallest region (in number of data points) that will be returned as a DMR. Default=4.
#' @param plotresult will output a graphical representation of each DMR to a pdf file in the specified directory.   Must specify build.
#' @return a list containing matrix of DMR results, as well as other values and parameters which can be passed to plotSatScanresult.  
#' @seealso \code{\link{SatScanrun}} which this function wraps, and \code{\link{plotSatScanresult}} which plots the results of this function. \code{\link{plotSatScanresult}} can be called within this function with plotresult=TRUE. 
#' @export
#' @examples
#' result<-methdiffSatScan(meth=meth, mydir=getwd(), xvalues="Index", splitchr=TRUE, splitcentromere=TRUE, build="hg18", smallestregion=4, plotresult=TRUE)
#' 
#' 




methdiffSatScan <- function(meth, mydir=NULL, xvalues=c("Index", "Position"), splitchr=FALSE, splitcentromere=FALSE, build=c("hg18", "hg19", "hg38"), smallestregion=4, plotresult=FALSE, ... ){
  
  library("rsatscan")
  library(methylSig)
  
  ##################### get the data in the right format
  
  if(class(meth)=="methylSigData"){
    cat("Data is a methylSigData object.")
  }else{
    stop("Data is not in a methylSigData object.  Please see the methylSigReadData function in the methylSig package.  ")
  }
  
  ##################### get/make the directory 
  if(is.null(mydir)){mydir<-getwd()}
  dir.create(paste(mydir, "/methdiffSatScan_tmp", sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  tmpdir<-paste(mydir, "/methdiffSatScan_tmp", sep="")
  
  
  ##################### 1: splitting by chromosome and centromere.
  #####################if we want to split by chromosome, see which chromosomes are present in the data
  if( splitchr==TRUE & splitcentromere==TRUE){
    
    
    ##################### if we want to split at centromere, read in centromeres for proper build
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
    if (length(badchrs)>0){
      stop(paste("Cannot split at the centromere for chromosome(s):\n",badchrs, ". \n  Please specify a chromosme as chr1, etc, or do not split by centromere."))
      
    }
    
    if (length(mychrlistorder)==length(mychrindlist)){
      mychrlistfinal<-mychrlistorder
    } else {
      mychrlistfinal<-mychrindlist
    }
    
    for (i in 1:length(mychrlistfinal)){
      mychr<-mychrlistfinal[i]
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
    
    
    
  } # end centromere and chr
  
  
  ##################### 2: splitting by chromosome, not centromere.
  
  if( splitchr==TRUE & splitcentromere==FALSE){
    found<-vector()
    loglik<-vector()
    
    mychrindlist<-rownames(table(meth@data.chr)) # this is not in order (1, 10, 11,...) instead of (1,2,..)
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
  
  ##################### 3: not splitting by chromosome or centromere
  if( splitchr==FALSE & splitcentromere==FALSE){
    
    result<-SatScanrun(meth=meth_chr, xvalues=xvalues, smallestregion=smallestregion, dir=tmpdir)
    found<-result$resultmat
    loglik<-result$logliks
    
      
  } # end no chr, no centromere
  
  ##################### 4: bad combination
  if( splitchr==FALSE & splitcentromere==TRUE){
    stop("Must split by chromosome to split by centromere.")
  } # end no chr, yes centromere
  
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







#' Plot the results of methdiffSatScan
#'
#' @param result is the result of running methdiffSatScan- a list containing a vector of likelihood values for each point, the DMRs which were found with p-value <0.05, the build, and the directory to which we want to write a pdf.
#' @param plottopn allows for the plotting of only the top n results by p-value, in case of a large number of identified DMRs.
#' @param plotpval allows for plotting of only DMRs with a p-value less than a given value (e.g. 0.001), in case of a large nunmber of identified DMRs.
#' @return nothing; results are output to pdf.  
#' @seealso \code{\link{methdiffSatScan}} which wraps this function 
#' @export
#' @examples
#' plotSatScanresult(result)
#' 
#' 
#' 


plotSatScanresult<-function(result, plottopn=NULL, plotpval=NULL){
  
  
  ########## this part controls how many results we want to plot.
  if (is.null(plottopn) & is.null(plotpval)){
    cat("plotting results to pdf")
    DMRs<-result$DMRs
  }
  if (!is.null(plottopn) & is.null(plotpval)){
    # order DMRs by pvalue.
    DMRs.pvalorder<-result$DMRs[with(result$DMRs, order(pval)), ]
    # get the p-val for the number of results specified
    mypval<-DMRs.pvalorder$pval[plottopn]
    # plot all with p-vals below this (for ties!)
    DMRs<-DMRs.pvalorder[DMRs.pvalorder$pval<=mypval,]
    
    if (dim(DMRs)[1]==plottopn){
      cat("plotting top", plottopn ,"results to pdf")
    }
    
    if (dim(DMRs)[1]>plottopn){
      cat("plotting top", dim(DMRs)[1] ,"results to pdf, to include those with equivalent p-values as top", plottopn)
    }
    
  }
  if (is.null(plottopn) & !is.null(plotpval)){
    cat("plotting results with p-value <=", plotpval ,"to pdf")
    DMRs<-result$DMRs[result$DMRs$pval<=plotpval,]
  }
  if (!is.null(plottopn) & !is.null(plotpval)){
    stop("Please set plottopn and/or plotpval to NULL.")
  }
  
  
  ########## resume with the function
  likelihoods<-result$likelihoods
  build<-result$build
  dir<-result$dir
  
  
  library(rtracklayer)
  library("GenomicFeatures")
  library("Gviz")
  
  
  downloadgenes<-makeTxDbFromUCSC( # this is a TxDb object, which is easy to plot with Gviz.
    genome=build,
    tablename="knownGene",
    transcript_ids=NULL,
    circ_seqs=DEFAULT_CIRC_SEQS,
    url="http://genome.ucsc.edu/cgi-bin/",
    goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
    miRBaseBuild=NA
  )
  
  
  ######## promoters from genes
  my.promoters<-trim(suppressWarnings(promoters(downloadgenes, upstream=2000, downstream=200)), use.names=TRUE)
  
  DMRs$chrnum<-as.numeric(gsub("chr", "", DMRs$chr))
  RSatScanDMRs.order<-with(DMRs, DMRs[order(chrnum, pos.start),])
  
  
  
  mychrindlist<-rownames(table(RSatScanDMRs.order$chr)) # this is not in order (1, 10, 11,...) instead of (1,2,..)
  mychrlistorder<-intersect(c(paste("chr",1:22, sep=""),"chrX", "chrY"), mychrindlist )
  
  pdf(file = paste(dir, "/RSatSscan_testplot.pdf", sep=""))
  
  DMRcount<-0
  for (i in 1:length(mychrlistorder)){
    mychr<-mychrlistorder[i]
    chrDMRs<-which(RSatScanDMRs.order$chr==mychr)
    
    
    ##########################    + ideogram     ##########################  
    
    itrack <- IdeogramTrack(genome = build, chromosome = mychr) # this makes the ideogram track, given build and chromosome
    
    ##########################    + genes     ##########################  
    
    grtrack <- GeneRegionTrack(downloadgenes, genome = build,
                               name = "Genes", showId=FALSE, #geneSymbols = TRUE, 
                               chromosome=mychr
                               #from=region$start-5000, to=region$end+5000,
                               #transcriptAnnotation = "transcript" 
                               #,collapseTranscripts="meta"
                               )
    
    #TranscriptDB does not have symbols?
    library("org.Hs.eg.db")
    symbols <- suppressWarnings(unlist(mapIds(org.Hs.eg.db, gene(grtrack), "SYMBOL", "ENTREZID", multiVals = "first")))
    
    
    
    
    a<-symbol(grtrack)
    b<-gene(grtrack)
    c<-symbols[gene(grtrack)]
    
    symbol(grtrack) <- symbols[gene(grtrack)]
    symbol(grtrack)[is.na(symbol(grtrack))]<-a[is.na(symbol(grtrack))]
    gene(grtrack)<- symbol(grtrack)  # this is new!
    ##########################    + promoters     ##########################  
    
    my.promoters.chr<-my.promoters[seqnames(my.promoters) == mychr  ]
    ptrack <- AnnotationTrack(my.promoters.chr, name = "Promoters")
    
    
    
    ##########################    + CpG islands     ##########################  
    
    CpGtrack=UcscTrack(track="CpG Islands", genome=build, chromosome=mychr, 
                       start="chromStart", end="chromEnd", name="CpG Islands", fill="#ADE6BA")
    
    
    
    for (DMRind in chrDMRs){  #
      RSatScanDMR<-RSatScanDMRs.order[DMRind,]
      DMRcount<-DMRcount+1
      cat("\n\n\n\n\n\nDMR", DMRcount, "of", dim(RSatScanDMRs.order)[1], "...\n\n\n\n\n\n")
      
   
      ##########################    + region     ##########################  
      
      region<-data.frame(RSatScanDMR["pos.start"],RSatScanDMR["pos.stop"], mychr)
      plotborder<-round((RSatScanDMR["pos.stop"]-RSatScanDMR["pos.start"])*1.0)
      plotstart<-max((RSatScanDMR["pos.start"]-plotborder), 0)
      plotstop<-min((RSatScanDMR["pos.stop"]+plotborder), end(tail(itrack@range@ranges, 1))  )
      
      
      names(region)<-c("start", "end", "chromosome")
      atrack <- AnnotationTrack(region, name = "RSatScan",  fill="#9966ff")
      
      
      ##########################    + likelihoods     ##########################  
      # this one is if we only want logliks in the region
      #mylogliks<-likelihoods[which(as.character(likelihoods$chr)==as.character(region$chromosome) & likelihoods$pos>=region$start  & likelihoods$pos<=region$end),]
      
      # this is if we want all logliks in the plotting area.
      mylogliks<-likelihoods[which(as.character(likelihoods$chr)==as.character(region$chromosome) & likelihoods$pos>=plotstart  & likelihoods$pos<=plotstop),]
      
      dtrack <- DataTrack(data = mylogliks$loglik, start = mylogliks$pos, end = mylogliks$pos, 
                          chromosome = as.character(region$chromosome), genome = build,
                          name = "LogLik", fill="#808080")

      
      ##########################    + gene names     ##########################  
      #grtrack@range@ranges #72163
      #symbol(grtrack)
      mygeneinds<-which(grtrack@range@ranges@start<=plotstop & (grtrack@range@ranges@start + grtrack@range@ranges@width -1 )>=plotstart)
      mygenes<-symbol(grtrack)[mygeneinds]
      namestarts<-pmax(plotstart+((plotstop-plotstart)*.12), (grtrack@range@ranges@start[mygeneinds]))
      
      mygenes.unique<-mygenes[!duplicated(mygenes)]
      namestarts.unique<-namestarts[!duplicated(mygenes)]
      
      ucgeneind<-sum(is.na(gene(grtrack)[mygeneinds]))
      
      if (length(namestarts)>0){
      nameTrack <- AnnotationTrack(start = namestarts.unique,
                                width = ((plotstop-plotstart)*.1), chromosome = as.character(region$chromosome), group = mygenes.unique,
                                genome = build, name = "Gene Names", showId = TRUE, 
                                col=NULL, 
                                fill="white")
      } else {nameTrack<- AnnotationTrack(start = plotstart,
                              width = ((plotstop-plotstart)*.1), chromosome = as.character(region$chromosome), id = "empty",
                              genome = build, name = "Gene Names", showFeatureId = FALSE, 
                              col=NULL, 
                              fill="white") 
      ucgeneind<-1
      }

      
      ##########################    plot all     ##########################  
      #failtracks<-list(itrack,  atrack, dtrack, grtrack, ptrack, CpGtrack )
      #succeedtracks<-list(itrack,  atrack, dtrack, grtrack, ptrack, CpGtrack )
      if(ucgeneind==0){
      plotTracks(list(itrack,  atrack, dtrack, grtrack, nameTrack, ptrack, CpGtrack ), from=plotstart, to=plotstop, geneSymbols = FALSE,
                          fontsize=15, type = "histogram", 
                         collapseTranscripts="meta" , 
                          main=paste(RSatScanDMR$chr, ":",RSatScanDMR$pos.start, "-",RSatScanDMR$pos.stop ), 
                          cex.main=0.75
      )
      }
      if(ucgeneind>0){
        plotTracks(list(itrack,  atrack, dtrack, grtrack, nameTrack,  ptrack, CpGtrack ), from=plotstart, to=plotstop, geneSymbols = FALSE,
                   fontsize=15, type = "histogram", 
                  main=paste(RSatScanDMR$chr, ":",RSatScanDMR$pos.start, "-",RSatScanDMR$pos.stop), 
                   cex.main=0.75  
                   #, showId = TRUE
                   #, add=TRUE 
        ) }
      
      
      
    } # end DMR in chr
  } # end chr
  
  dev.off()
  
} # end function

#save(failtracks, file="H:/Methylation/Methods_comparison/methdiffSatScan_testing_code/failtracks.RData")
#save(succeedtracks, file="H:/Methylation/Methods_comparison/methdiffSatScan_testing_code/succeedtracks.RData")

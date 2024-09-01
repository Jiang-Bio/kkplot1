

#' plot_track
#'
#' @param Grange object
#' @param show_peak Grange object
plot_track <- function(zoom_region,bwfile,y_max=40,track_color="#80B1D3",peak_color="#FB8072",sitename="TF",
                       show_gene=TRUE,show_peak=NULL,genome=NULL,ylables=c("10","40"),tick.pos=c(10,40),
                       label_margin=0.08,show_tatget_gene=NULL){
  if (is.null(genome)) {
    genome <- readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/ara.fa")}

  genome <- GRanges(seqnames = names(genome),
                    ranges = IRanges(start = rep(1, length(genome)), width = width(genome)),
                    strand = rep("*", length(genome)))

  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  if (show_gene) {
    txdbbb <- TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51
    seqlevels(txdbbb, pruning.mode="coarse") <- seqlevels(BSgenome.Athaliana.TAIR.TAIR9::Athaliana)
    kp <- plotKaryotype(zoom = zoom_region,genome = genome,plot.type =4,plot.params = pp)
    genes.data <- makeGenesDataFromTxDb(txdbbb,
                                        karyoplot=kp,
                                        plot.transcripts = TRUE,
                                        plot.transcripts.structure = TRUE)
    genes.data <- mergeTranscripts(genes.data)
    if (!is.null(show_tatget_gene)) {
      message(genes.data$genes$gene_id)
      num <- which(genes.data$genes$gene_id==show_tatget_gene)
      genes.data$genes <- genes.data$genes[num]
      genes.data$transcripts <- genes.data$transcripts[num]
      genes.data$coding.exons <- genes.data$coding.exons[num]
      genes.data$non.coding.exons <- genes.data$non.coding.exons[num]
    }

    kpPlotGenes(kp, data=genes.data, r0=0, r1=0.05,add.strand.marks = TRUE,mark.height = 0.5,
                introns.col = "gray",mark.width = 1,mark.distance = 0.1,marks.col = "gray")
  }else {
    kp <- plotKaryotype(zoom = g,genome = genome,plot.type =4,plot.params = pp)
  }
  kpAddBaseNumbers(kp, tick.dist = 1000,tick.len = 2.5 ,minor.tick.len = 1,minor.tick.dist = 200,
                   add.units = TRUE,units = "kb" ,cex=0.5, digits = 2)

  if (!is.null(show_peak)) {
    peakinfo <-  show_peak
    kpPlotRegions(kp, peakinfo, col=peak_color, r0=0.15, r1=0.18)
    kpAddLabels(kp, labels = glue("{sitename} \n binding site"), r0=0.15, r1=0.18, cex = 0.8)
  }

  for(i in seq_len(length(bwfile))) {
    bigwig.file <-bwfile[i]
    at <- autotrack(i, length(bwfile), r0=0.2, r1=1,margin = 0.1)
    kp <- kpPlotBigWig(kp, data=bigwig.file, ymax=y_max,r0=at$r0, r1=at$r1,col = track_color)
    computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    kpAxis(kp, r0=at$r0, r1=at$r1,cex=0.5, ymin=0, ymax=y_max,tick.pos = tick.pos,
           labels = ylables )
    kpAddLabels(kp, labels = names(bwfile)[i], r0=at$r0, r1=at$r1,cex=0.5,srt=0,
                pos=3, offset = 10,label.margin = label_margin)
  }
}

#' Title
#'
#' @param peak narrowPeak Grange object
#' @param top_peak_num num
#' @param go_num How many are displayed
peak_to_GO <- function(peak,top_peak_num=2000,tssregion=c(-1000,1000),go_num=30){
  df <- import(i)
  df <- df[grep("^Chr", df@seqnames),]
  df <- df[order(df$score,decreasing = TRUE)][1:top_peak_num]
  seqlevels(df, pruning.mode="coarse") <- seqlevels(TxDb.Athaliana.BioMart.plantsmart51)
  df <- ChIPseeker::annotatePeak(peak = df,TxDb = TxDb.Athaliana.BioMart.plantsmart51::TxDb.Athaliana.BioMart.plantsmart51,
                                 tssRegion = tssregion,annoDb =org.At.tair.db) %>% as.data.frame()
  genes <- df$geneId[grep("Promoter",df$annotation)] %>% unique()
  ego <- enrichGO(genes, OrgDb = org.At.tair.db::org.At.tair.db, keyType = "TAIR", ont = "ALL", pvalueCutoff =1,qvalueCutoff = 1) %>%
    as.data.frame()
  top=go_num
  Description <- ego[1:top,]$Description %>% unique()
  df <- data.frame("Description"=Description) %>% left_join(ego [,c("Description","pvalue","Count","GeneRatio")]) %>%
    set_colnames(c("Description","pvalue","Count","GeneRatio"))
  df$pvalue <-  -log10(df$pvalue)
  df$Description <- factor(df$Description, levels = rev(df$Description)) #设置不按照首字母大小排序
  df$Description <- rev(df$Description) ;df$pvalue <- rev(df$pvalue) ;df$Count <- rev(df$Count)
  plot <- ggplot(df, aes(y = Description    ,x = Count , fill = pvalue )) +
    geom_bar(stat = "identity") +
    scale_fill_gradientn(colours =brewer.pal(5,"Blues"),name="-log10(pvalue)" ) +
    labs(title = "", x = "Category", y = "Value") +
    theme_bw()+
    geom_text(aes(x = 0,label =Description,hjust ="left"))+
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 15),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          legend.text = element_text(size = 10),
          legend.title  = element_text(size = 10),
          legend.position = "right")+  xlab("Count")+
    ylab("")
  return(plot)

}

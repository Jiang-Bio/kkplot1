#' Title
#'
#' @param pwmlist
#' @param outfile
#' @param weigh
#'
#' @return
#' @export
#'
#' @examples
TF_logo <- function(pwmlist,outfile,weigh=NULL){
  motif_pics=BiocParallel::bplapply(1:length(pwmlist),function(x){
    t_pwmlist <- universalmotif::convert_motifs(pwmlist[[x]])
    motif_pic=suppressMessages(universalmotif::view_motifs(t_pwmlist,show.positions = F,tryRC = T,
                                                           min.mean.ic = 0,fit.to.height = 1,names.pos = "right",
                                                           normalise.scores = TRUE,min.position.ic = 0,text.size = 4,min.height = 0,min.overlap = 20)+
                                 ylab("")+
                                 scale_y_continuous() +
                                 theme(axis.ticks.y =element_blank(),
                                       axis.line.y = element_blank(),
                                       axis.text.y = element_blank())  )
    ggplotGrob(motif_pic)
  },BPPARAM = BiocParallel::MulticoreParam(workers = 30))

  if (length(pwmlist)>100) {
    blank_plot <- ggplot() +
      theme_void()+
      xlim(0, 6) +
      ylim(0, 100)#)#+xlim(0,2)
    anno_plots <- c()
    anno_texts <- c()
    for (i in 1:100) {
      cmd <- glue("annotation_custom(motif_pics[[{i}]],xmin =0.5,xmax={0.5+
                ncol(pwmlist[[i]]@profileMatrix)*0.1},ymin={99-i},ymax={99-(i-2.5)})")
      cmd2 <- glue("annotate('text', x = 0.4, y = {100.5-(i)}, label = '{pwmlist[[i]]@name}',size = 2)")
      anno_plots <- append(anno_plots,cmd)
      anno_texts <- append(anno_texts,cmd2)
    }
    for (i in 101:length(pwmlist)) {
      cmd <- glue("annotation_custom(motif_pics[[{i}]],xmin =3.5,xmax={3.5+
                ncol(pwmlist[[i]]@profileMatrix)*0.1},ymin={99-i+100},ymax={99-(i-2.5)+100})")
      cmd2 <- glue("annotate('text', x = 3.4, y = {100.5-(i)+100}, label = '{pwmlist[[i]]@name}',size = 2)")
      anno_plots <- append(anno_plots,cmd)
      anno_texts <- append(anno_texts,cmd2)
    }
    anno_plots %<>% paste0(collapse = "+")
    anno_texts %<>% paste0(collapse = "+")
    plot_text=paste0("blank_plot+ ",anno_plots,"+",anno_texts)
    plot_1 <- eval(parse(text = plot_text))
    ggsave(filename = outfile,plot = plot_1,width = 5,height = 20)
  }else{
    if (is.null(weigh)) {
      blank_plot <- ggplot() +
        theme_void()+
        xlim(0, 3) +
        ylim(0, 100)
    }else{
      blank_plot <- ggplot() +
        theme_void()+
        xlim(0, weigh) +
        ylim(0, 100)
    }

    anno_plots <- c()
    anno_texts <- c()
    for (i in 1:length(pwmlist)) {
      cmd <- glue("annotation_custom(motif_pics[[{i}]],xmin =0.5,xmax={0.5+
                length(pwmlist[[i]]@profileMatrix[1,])*0.1},ymin={99-i},ymax={99-(i-2.5)})")
      cmd2 <- glue("annotate('text', x = 0.4, y = {100.5-(i)}, label = '{pwmlist[[i]]@name}',size = 1)")
      anno_plots <- append(anno_plots,cmd)
      anno_texts <- append(anno_texts,cmd2)
    }
    anno_plots %<>% paste0(collapse = "+")
    anno_texts %<>% paste0(collapse = "+")
    plot_text=paste0("blank_plot+ ",anno_plots,"+",anno_texts)
    plot_1 <- eval(parse(text = plot_text))
    ggsave(filename = outfile,plot = plot_1,width = 5,height = 20)
  }

}

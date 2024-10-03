#' Title
#'
#' @param SRR
#' @param outpath
#' @param Private_key
#'
#' @return
#' @export
#'
#' @examples
download_data_ascp <-function(SRR,outpath,Private_key="/wrk/jiangdingkun/.aspera/connect/etc/asperaweb_id_dsa.openssh"){
  mclapply(1:length(SRR),function(i){
    srrxx <- SRR[i] %>% substr(start = 1, stop = 6)
    # RUN LibraryLayout
    # 1 SRR7663611        PAIRED
    # 2 SRR7663612        PAIRED
    str <- strsplit(SRR[i], "")[[1]]
    if (length(str)==11) {
      id <- str[10:11] %>% paste0(collapse = "") %>% glue("0",.)
    }else{
      id <- tail(str, 1) %>% glue("00",.)
    }

    if (!dir.exists(oupath)) {
      dir.create(glue("{oupath}/SINGLE"))
      dir.create(glue("{oupath}/PAIRED"))}

    if (SRR$LibraryLayout[i]=="SINGLE") {
      cmd <- glue("ascp -i {Private_key} -l 200M -P 33001 -vQT -k 1 era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srrxx}/{id}/{SRR[i]}/{SRR[i]}.fastq.gz {oupath}/SINGLE/.")
      if (file.exists(glue("{oupath}/{SRR[i]}.fastq.gz.aspx"))&&file.exists(glue("{oupath}/{SRR[i]}.fastq.gz"))) {
        system(cmd,intern = T)
      }
      if (!file.exists(glue("{oupath}/SINGLE/{SRR[i]}.fastq.gz"))) {
        system(cmd,intern = T)
      }
    }
    if (SRR$LibraryLayout[i]=="PAIRED") {
      message(glue("the file {SRR[i]} is downloading"))
      cmd <- glue("ascp -i {Private_key} -l 200M -P 33001 -vQT -k 1 era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srrxx}/{id}/{SRR[i]}/{SRR[i]}_1.fastq.gz {oupath}/PAIRED/.")
      if (file.exists(glue("{oupath}/PAIRED/{SRR[i]}_1.fastq.gz.aspx"))&&file.exists(glue("{oupath}/PAIRED/{SRR[i]}_1.fastq.gz"))) {
        system(cmd,intern = T)
      }
      if (!file.exists(glue("{oupath}/PAIRED/{SRR[i]}_1.fastq.gz"))) {
        system(cmd,intern = T)
      }

      cmd <- glue("ascp -i {Private_key} -l 200M -P 33001 -vQT -k 1 era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/{srrxx}/{id}/{SRR[i]}/{SRR[i]}_2.fastq.gz {oupath}/PAIRED/.")
      if (file.exists(glue("{oupath}/PAIRED/{SRR[i]}_2.fastq.gz.aspx"))&&file.exists(glue("{oupath}/PAIRED/{SRR[i]}_2.fastq.gz"))) {
        system(cmd,intern = T)
      }
      if (!file.exists(glue("{oupath}/PAIRED/{SRR[i]}_2.fastq.gz"))) {
        system(cmd,intern = T)
      }
    }
  },mc.cores = 10)
}




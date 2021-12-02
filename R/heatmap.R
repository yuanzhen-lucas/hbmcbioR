#'you can make a heatmap to see your fasta file kmer
#'
#'@param sampleinfor a tibble or datafram include the information of batch well tissue ligand length and file dress
#'@param outpath to save the results
#'@return different ligand length heatmap png and pdf
#'
#'
#'@example
#' sample <- readxl::read_excel("/wrk/yuanzhen/heatmap/202109_root_ATIs.xlsx")
#' heatmap_kmer(sample,"/wrk/yuanzhen/heatmap")
#'




heatmap_kmer <- function(sampleinfor,outpath){
  kmer <- lapply(sample$file,function(x){
    fjComm::kmerCntBit(strings =fjComm::getSeq_fqfachrFile(x), k = 9L, diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = 10)
  })
  saveRDS(kmer,paste0("kmerfile",".Rds"))

  kmer_cn <- lapply(kmer, function(x){
    arrange(x,desc(counts)) %>% filter(row_number() <= 1000) %>% mutate(log_counts = log(counts)) %>% select(-(counts))
  })

  kmer_merge <- kmer_cn %>% purrr::reduce(full_join, by = "kmer")

  kmer_merge[is.na(kmer_merge)] <- 1

  name <- paste0(sample$well,"_",sample$`ligand length`,"_",sample$batch)

  m <- as.matrix(kmer_merge[, -1])
  rownames(m) <- kmer_merge$kmer
  colnames(m) <- name
  mdt <- as.data.frame(m)
  m30 <- dplyr::select(mdt,contains("30"))
  m101 <- dplyr::select(mdt,contains("101"))

  png30 <- pheatmap::pheatmap(m30 %>% as.matrix(),show_rownames = F)
  png101 <- pheatmap::pheatmap(m101 %>% as.matrix(),show_rownames = F)

  ggplot2::ggsave(png30,paste0(outpath,"/","30n_kmer",".png"))
  ggplot2::ggsave(png30,paste0(outpath,"/","30n_kmer"),device="pdf")
  ggplot2::ggsave(png101,paste0(outpath,"/","101n_kmer",".png"))
  ggplot2::ggsave(png101,paste0(outpath,"/","101n_kmer",".pdf"))

}


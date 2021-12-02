#' ATACseq peak annotation
#' @param peaks the peak call file from the fseq2 of GRanges object
#' @param txdb  you interested species annotation genome file of TxDb object
#' @param tss_region the number of base pair from tss ,default is 3000
#' @param prefix the name for save your result
#' @param tfname the TF family name can be plot with RNA-seq and ATAC-seq,default is WRKY
#'
#'
#' @example
#' genomefile <- readDNAStringSet("/wrk/data/genome/yz_genome_data/aragenome/TAIR10_chr_all.fas")
#' txdb <- makeTxDbFromGFF("/wrk/data/genome/yz_genome_data/aragenome/Athaliana.gff3")
#' peaks <- ChIPseeker::readPeakFile("/wrk/yuanzhen/atacseq/fseq2_result_peaks.narrowPeak")
#' rnaseqtxt <- "/wrk/yuanzhen/gene_tpm.txt"
#' atac_peakanno(peaks,txdb,rnaseqtxt,prefix="open_chr")


#library(GenomicRanges)
#library(ChIPseeker)
#library(GenomicFeatures)
#library(rGADEM)
#library(IRanges)
#library(magrittr)
#library(ggpubr)
#library(BSgenome.Athaliana10.TAIR.TAIR10)
#library(org.At.tair.db)
atac_peakanno <- function(peaks,txdb,rnaseqtxt,tss_region=3000,prefix,tfname="WRKY"){
  #peak annotation
  genes <- genes(txdb)
  if(all(seqlevels(peaks) %in% c("Chr1","Chr2","Chr3","Chr4","Chr5","chloroplast"))){
    seqlevels(peaks) <- c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrC")
  }else if(all(seqlevels(peaks) %in% c("Chr1","Chr2","Chr3","Chr4","Chr5","mitochondria"))){
    seqlevels(peaks) <- c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrM")
  }else{
    seqlevels(peaks) <- c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrC","ChrM")
  }
  colnames(peaks@elementMetadata) <- c("appname","score","strand","signalValue","pValue(-log10)","qValue(-log10)","peak")
  peaks@elementMetadata@listData$appname <- NULL
  peakAnno <- ChIPseeker::annotatePeak(peaks, tssRegion=c(-tss_region, tss_region), TxDb=txdb)
  write_rds(peakAnno@anno,paste0(prefix,"_","peak_annotation.RDS"))

#the coverage of bam file
  cov_raw <- coverage("/wrk/yuanzhen/atacseq/pre_propressed/trim_map/dmso10u_30min.sorted_rmdup.bam")
  vis_width=3000*2+1
  tss=promoters(genes,upstream = 0,downstream = 1)
  tss=resize(tss,1000,fix = "end") %>% trim()
  tss=tss[width(tss)==1000]
  rngs=ranges(tss)
  start_=IRanges::start(rngs)
  filter_=start_>0
  tss=tss[filter_]
  seqlevels(tss) %<>% str_replace("Chr","chr")
  seqlevels(tss,pruning.mode="coarse")=seqlevels(tss)[seqlevels(tss) %>% str_detect("[0-9]$")]
  names(cov_raw) %<>% str_replace("Chr","chr")
  for (i in 1:length(cov_raw)){
    cov_raw[[i]]=c(cov_raw[[i]],Rle(0L,10000))
  }
  cov <- cov_raw[tss]
  cov %<>% unlist(use.names = F) %>% as.integer()
  cov_mt=matrix(cov,ncol=1000,byrow = T)
  cov_sum=data.frame(rowSums(cov_mt))
  colnames(cov_sum) <- "coverageSum"
  cov_sum$gene_id <- tss@elementMetadata$gene_id

#the Correlation of RNA-seq and ATAC-seq,you can choose a specific TF family
  gene_tpm <- read.delim(rnaseqtxt)
  gene_exp <- data.frame(dmso_exp=(gene_tpm$dmso1 %>% as.numeric() + gene_tpm$dmso2 %>% as.numeric()) / 2)
  gene_exp$gene_id <- rownames(gene_tpm)
  gen_total <- cov_sum %>% inner_join(gene_exp,by="gene_id")
  gen_total$dmso_exp[which(gen_total$dmso_exp < 1)] <- 1
  gen_total$coverageSum[which(gen_total$coverageSum < 1)] <- 1
  gen_total$log_cov <- log(gen_total$coverageSum)
  gen_total$log_gexp <- log(gen_total$dmso_exp)
  gen_filter <- filter(gen_total,log_cov > 5 & log_gexp >0)

  #the TF family gene id from org.At.tair.db database
  dt <- AnnotationDbi::as.data.frame(AnnotationDbi::select(org.At.tair.db, keys=keys,columns="SYMBOL", keytypes="GENENAME"))
  x_tn <- paste0("^",tfname)
  dt1=dt %>% filter(!is.na(SYMBOL))
  ft <- grepl(pattern = x_tn,dt1$SYMBOL)
  dt_wk <- dt1$SYMBOL[ft] %>% as.data.frame()
  colnames(dt_wk) <- "SYMBOL"
  dt_name <- left_join(dt_wk,dt1,by="SYMBOL")
  colnames(dt_name) <- c("family","gene_id")
  fi_wk <- left_join(dt_name,gen_filter,by="gene_id",rm.) %>% drop_na()
  p2 <- ggplot(gen_filter, aes(x=log_cov, y=log_gexp)) +
    geom_point(size=0.1)+
    geom_smooth(method = 'lm', se = F, color = 'red')+
    theme_bw()+
    stat_cor(data=gen_filter, method = "spearman")+
    scale_x_continuous(limits = c(8, 11))+
    scale_y_continuous(limits = c(1, 10))
  ggsave(p2)


  #de novo motif with rGADEM
  rang <- ranges(peaks)
  bed <- data.frame(chr=as.factor(seqnames(peaks)),start=as.numeric(rang@start),end=as.numeric(rang@start + rang@width -1))
  seqbio <- Biostrings::DNAStringSet()
  for(a in 1:length(levels(seqnames(peaks)))){
    chr <- bed %>% filter(chr==levels(seqnames(peaks))[[a]])
    seqbio <- c(seqbio,Biostrings::DNAStringSet(BSgenome.Athaliana10.TAIR.TAIR10[[levels(seqnames(peaks))[[a]]]], start=chr$start, end=chr$end))
  }
  gadem <- GADEM(seqbio,verbose=1,genome=BSgenome.Athaliana10.TAIR.TAIR10)
  write_rds(gadem,paste0(prefix,"_","motif.RDS"))
}





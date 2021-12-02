    #' just compare max number
    #'
    #' @param v
    #'
    #' @return
    #' @export
    #'
    #' @examples
    .getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }



#' motif count number in a seqs
#' @description put a DNa sequence file,then you can get ecah motif count number
#' @param seqs a DNA sequence file(just include base)
#' @param pwmlist motifs of PWMartixList format
#'
#' @return
#' @export
#'
#' @examples
#' motif_count_ <- motif_counts_(seqs,pwmlist)
motif_counts_ <- function(seqs,pwmlist){
  dnabase <- readLines(seqs)
  dnaseq <- DNAStringSet(dnabase)
  motif_ix_SummarizedExperiment <- matchMotifs(pwmlist,dnabase,bg="even",out = "scores")
  motif_ct <- motifmatchr::motifCounts(motif_ix_SummarizedExperiment) %>% Matrix::colSums()
  motif_ffr <- motif_ct / length(dnaseq)

  test <- Biostrings::alphabetFrequency(dnaseq, baseOnly=TRUE, collapse=TRUE)
  seqnons <- create_sequences(seqnum = length(dnaseq), seqlen = .getmode(width(dnaseq)),
                              freqs = c(A=test[[1]] / sum(test), C=test[[2]] / sum(test),
                                        G=test[[3]] / sum(test), T=test[[4]] / sum(test)))
  motif_ix_SummarizedExperiment_nomr <- matchMotifs(pwmlist,seqnons,bg="even",out = "scores")
  motif_ct_nomr <- motifmatchr::motifCounts(motif_ix_SummarizedExperiment_nomr) %>% Matrix::colSums()
  motif_ffr_nomr <- motif_ct_nomr / length(seqnons)

  motif_deno <- motif_ffr - motif_ffr_nomr
  fil_motif <- motif_deno <= 0
  motif_deno[fil_motif] <- 0.000001
  motif_deno
}


#' match your motif to the most similar known motif form database
#' @description this return an TFBSTools format motif ,it's ID includes the name,family and similarty of known motif
#' @seealso  [motif_database_html()],to see which known motif include
#' @param PFMatrix
#'
#' @return
#' @export
#'
#' @examples
motif_get_name <- function(PFMatrix){
  motifdatabase <- read_motifs("/wrk/yuanzhen/motifdb/motif_database")
  pwm_library <- read_rds("/wrk/yuanzhen/motifdb/pwm_library.Rds")

  motiflogo <- universalmotif::convert_motifs(PFMatrix)
  # find the most similar motif to our motif
  pwm_sim <- PWMSimilarity(pwm_library, PFMatrix, method = 'Pearson')
  pwm_sim_rc <- PWMSimilarity(pwm_library, convert_motifs(motif_rc(motiflogo),"TFBSTools-PWMatrix"), method = 'Pearson')
  pwm_sim <- pmax(pwm_sim,pwm_sim_rc)

  # extract the motif names from the pwm library
  pwm_library_list = lapply(pwm_library, function(x){
    data.frame(ID = TFBSTools::ID(x))})

  # combine the list into one data frame
  pwm_library_dt = dplyr::bind_rows(pwm_library_list)

  # fetch the similarity of each motif to our unknown motif
  pwm_library_dt$similarity = pwm_sim

  # find the most similar motif in the library
  pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]
  pwm_library_dt$col = rownames(pwm_library_dt)
  x_n <- head(pwm_library_dt,1)
  ge <- readr::read_rds("/wrk/yuanzhen/motifdb/motif_gene.Rds")
  name_ <- lapply(motifdatabase,function(x)x["name"]) %>% unlist()
  pwmnub <- x_n$col %>% as.numeric()
  pwm_x <-  pwm_library[[pwmnub]]
  know_motif <- motifdatabase[[which(name_==pwm_x@ID)]]
  if(str_detect(pwm_x@ID,"_cisbp")){
    tf_locus <- pwm_x@tags$DBID
    database_source <- "CIS-BP"
  }
  if(str_detect(pwm_x@ID,"_plantDB")){
    tf_locus <- pwm_x@tags$Gene_id
    database_source <- "plantDB"
  }

  if(str_detect(pwm_x@ID,"_plantpan")){
    tf_locus <- pwm_x@tags$TF_Locus
    database_source <- "plantPAN"
  }

  name <- ge %>% dplyr::filter(gene==tf_locus)

  if(nrow(name)==1){
    TF_name <- name$TF_name
    TF_family <- name$TF_family
  }else if(nrow(name)>1){
    TF_name <- name$TF_name[[1]]
    TF_family <- name$TF_family[[1]]
  }else{
    TF_name <- "NO_source_name"
    TF_family <- "NO_source_family"
  }
  x_name_motif <- paste0(TF_name,"_",TF_family,"_",round(x_n$similarity,3))
  x_name_motif <- str_replace(x_name_motif,"/","&")
  PFMatrix@ID <- x_name_motif
  PFMatrix@name <- x_name_motif
  PFMatrix
}


#' attach motif pictures to plot
#'
#' @param pwmlist motif list
#' @param tdt30
#' @param tdt101
#' @param tdt30
#' @param yz_30
#' @param yz_101
#' @param threads
#' @return
#' @export
#' @seealso [info_filter()]
#' @examples
pic_attach <- function(pwmlist,tdt30,tdt101,yz_30,yz_101,threads){
  dir_temp <- getwd()
  dir_index <- paste0(dir_temp,"/motif_index")
  dir.create(dir_index)
  parallelParam= BiocParallel::MulticoreParam(workers = threads)
  bplapply(1:length(pwmlist),
           function(x){
             t_pwmlist <- convert_motifs(pwmlist[[x]])
             con_pwmlist <- convert_type(t_pwmlist,"ICM")
             width <- str_count(t_pwmlist["consensus"]) / 1.5
             height <- t_pwmlist@icscore / 5
             ggsave(paste0(dir_index,"/",t_pwmlist["name"],".png"),
                    view_motifs(t_pwmlist,show.positions = F)+ylab("")+theme(axis.ticks.y =element_blank(),axis.line.y = element_blank(),axis.text.y = element_blank())#labs(title = t_pwmlist["name"])+theme(plot.title = element_text(hjust = 0.5,vjust = -2))
                    ,width = width,height = height)}
           ,BPPARAM = parallelParam)

  file_mt <- Sys.glob(paste0(dir_index,"/","*.png"))
  labels_motif_30 <- c()
  for (i in 1:length(file_mt)){
    i_pwmlist <- convert_motifs(pwmlist[[i]])
    img.name <- i_pwmlist["name"]
    width <- str_count(i_pwmlist["consensus"]) / 1.5 * 4
    height <- i_pwmlist@icscore / 5 *4
    labels_motif_30 <- c(labels_motif_30, paste0("<img src='", file_mt[[i]],  "' width=","'",width,"'" ," ","height='",height,"'"," /><br>*", tdt30$tfname[[i]],"*"))
  }

  labels_motif_101 <- c()
  for (i in 1:length(file_mt)){
    i_pwmlist <- convert_motifs(pwmlist[[i]])
    img.name <- i_pwmlist["name"]
    width <- str_count(i_pwmlist["consensus"]) / 1.5 * 4
    height <- i_pwmlist@icscore / 5 *4
    labels_motif_101 <- c(labels_motif_101, paste0("<img src='", file_mt[[i]],  "' width=","'",width,"'" ," ","height='",height,"'"," /><br>*", tdt101$tfname[[i]],"*"))
  }

  plot_img_30 <- yz_30 + ggplot2::scale_x_discrete(name = NULL,
                                                   labels = rev(labels_motif_30)) +
    ggplot2::theme(axis.text.y = element_markdown(color = "black", size = 3,align_widths=T,align_heights=T))

  plot_img_101 <- yz_101 + ggplot2::scale_x_discrete(name = NULL,
                                                     labels = rev(labels_motif_101)) +
    ggplot2::theme(axis.text.y = element_markdown(color = "black", size = 3,align_widths=T,align_heights=T))

  ggsave("plot_30_motif.png",plot_img_30,device = "png",type ="cairo")
  ggsave("plot_101_motif.png",plot_img_101,device = "png",type ="cairo")
}





#' you can compare two fasta file include the enrichment of your interested motifs
#'
#'
#' @param sample_info a DNAStringSetList containing the reads must be fasta file
#' @param pwmlist  a PWMatrixList containing the pwm of motifs
#' @param by ligand or tissue
#' @param filtername your group of by
#' @param threads the number threads you use
#' @param
#' @return
#'a plot png and pdf file
#' @export
#'
#' @examples
#'
#' ### make pwmlist
#' motifs <- readRDS("/wrk/zhu/test/motif_curation/curated_30N_motif_20211014.rds")
#' pwmlist <- PWMatrixList()
#' for(f in 1:length(motifs)){
#' pwmlist[[f]] <- convert_motifs(motifs[[f]], class = "TFBSTools-PWMatrix")
#' }
#'
#' ##read sampleinfo
#' sample_info <- read_excel("/wrk/yuanzhen/all_ati/rati.xlsx")
#'
#' info_filter(sample_info,pwmlist,30,by="ligand",filtername=c("30n","101n"))
#'
info_filter <- function(sample_info,pwmlist,threads,by="tissue",filtername="root"){
  browser()
  parallelParam= BiocParallel::MulticoreParam(workers = threads)
  pwmlist_re <- bplapply(pwmlist,function(x){motif_get_name(x)},BPPARAM = parallelParam)
  pwmlist_new <- PWMatrixList()
  for(i in 1:length(pwmlist_re)){
    pwmlist_new[[i]] <- pwmlist_re[[i]]
  }

  tfname <- lapply(pwmlist_new,function(x)x@name) %>% unlist()
  parallelParam= BiocParallel::MulticoreParam(workers = threads)
  motif_ra <- bplapply(sample_info$file, function(x){motif_counts_(x,pwmlist_new)}, BPPARAM = parallelParam)
  motif_ra <- do.call(cbind,motif_ra)
  colnames(motif_ra) <- paste0(sample_info$batch,"__",sample_info$well,"__",sample_info$tissue,"_",sample_info$ligand_length)
  rownames(motif_ra) <- tfname


  if(by=="tissue"){
    dt_nor <- motif_ra %>% log() %>% as.data.frame()
    name_fit <- str_detect(colnames(dt_nor),filtername)
    name_dt <- dt_nor[name_fit]
    pheatmap::pheatmap(name_dt %>% as.matrix())
  }else if(by=="ligand"){
    dt_mot <- motif_ra  %>% as.data.frame()
    filter_30 <- str_detect(colnames(dt_mot),filtername[[1]])
    filter_101 <- str_detect(colnames(dt_mot),filtername[[2]])
    dt_30 <- dt_mot[filter_30]
    n101_dt <- dt_mot[filter_101]
    dt_101_3 <- n101_dt / 3

    #dt_30$sd=apply(dt_30[,1:ncol(dt_30)],1,sd)
    dt_30$mean=apply(dt_30[,1:ncol(dt_30)],1,mean)
    dt_30$filename <- rep("ATI",nrow(dt_30))
    dt_30$tfname <- rownames(dt_30)
    #dt_101_3$sd=apply(dt_101_3[,1:ncol(dt_101_3)],1,sd)
    dt_101_3$mean=apply(dt_101_3[,1:ncol(dt_101_3)],1,mean)
    dt_101_3$filename <- rep("ApTI",nrow(dt_101_3))
    dt_101_3$tfname <- rownames(dt_101_3)

    tdt1 <- dt_30 %>% dplyr::select(tfname,filename,mean)
    tdt2 <- dt_101_3 %>% dplyr::select(tfname,filename,mean)


    plot_ati <- ggplot2::ggplot(data = tdt1, mapping = aes(x=tfname,y = mean)) +
      coord_flip()+
      scale_y_continuous(expand=c(0,0))+
      fjComm::gg_theme_Publication() +
      geom_bar(stat = 'identity',width = 0.5,fill="red") +
      #geom_point(aes(x=filename,y=sample_value),color = sample)+
      #geom_errorbar( aes(x=tfname, ymin=mean-sd, ymax=mean+sd), width=0.9, colour="black", alpha=0.5, size=0.1)+
      #theme(strip.text.y = element_blank(),axis.text.y= element_blank(),axis.ticks = element_blank())+
      ylab("Motif frequency")+xlab("")+theme(axis.ticks = element_blank())

    plot_apti <- ggplot2::ggplot(data = tdt2, mapping = aes(x=tfname,y = mean)) +
      coord_flip() +
      scale_y_continuous(expand=c(0,0))+
      fjComm::gg_theme_Publication() +
      geom_bar(stat = 'identity',width = 0.5,position = 'dodge',fill="black") +
      #geom_point(aes(x=filename,y=sample_value),color = sample)+
      #geom_errorbar( aes(x=filename, ymin=mean-sd, ymax=mean+sd), width=0.9, colour="black", alpha=0.5, size=0.1)+
      #theme(strip.text.y = element_blank(),axis.text.y= element_blank(),axis.ticks = element_blank())+
      ylab("Motif frequency")+xlab("")+theme(axis.ticks = element_blank())


    pic_attach(pwmlist = pwmlist_new,tdt30 = tdt1,tdt101 = tdt2,yz_30 = plot_ati,yz_101 = plot_apti,threads=threads)

    ggsave("ati.png",plot_ati,device = "png")
    ggsave("ati.pdf",plot_ati,device = "pdf")
    ggsave("apti.png",plot_apti,device = "png")
    ggsave("apti.pdf",plot_apti,device = "pdf")
  }


}






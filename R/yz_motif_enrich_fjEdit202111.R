#' match your motif to the most similar known motif form database
#' @description this return an TFBSTools format motif ,it's ID includes the name,family and similarty of known motif
#' @seealso  [motif_database_html()],to see which known motif include
#' @param PFMatrix
#'
#' @return
#' @export
#'
#' @examples
#' pwmlist="/wrk/zhu/test/motif_curation/curated_30N_motif_pwmlist_20211014.rds" %>% readRDS()
#' motif_get_name(pwmlist[[1]])
motif_get_name <- function(PFMatrix){
  motifdatabase <- universalmotif::read_motifs("/wrk/yuanzhen/motifdb/motif_database")
  pwm_library <- read_rds("/wrk/yuanzhen/motifdb/pwm_library.Rds")

  motiflogo <- universalmotif::convert_motifs(PFMatrix)
  # find the most similar motif to our motif
  pwm_sim <- TFBSTools::PWMSimilarity(pwm_library, PFMatrix, method = 'Pearson')
  pwm_sim_rc <- TFBSTools::PWMSimilarity(pwm_library, universalmotif::convert_motifs(universalmotif::motif_rc(motiflogo),"TFBSTools-PWMatrix"), method = 'Pearson')
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
  x_n <- head(pwm_library_dt,1) # get top n similar TFs, set to 1 here
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
  x_name_motif <- paste0(TF_family,"_",TF_name,"_",round(x_n$similarity,3))
  x_name_motif <- str_replace(x_name_motif,"/","&")
  PFMatrix@ID <- x_name_motif
  PFMatrix@name <- x_name_motif
  PFMatrix
}

#' motif enrichment analysis for a target library
#' @description count motif occurances in a target file (or character vector) of sequences, the frequency is normalized to a background library with the same 1-nt frequency but got randomized
#' @param seqs a DNA sequence file or a character vector of sequences
#' @param pwmlist motifs of PWMartixList format
#'
#' @return
#' @export
#'
#' @examples
#' pwmlist="/wrk/zhu/test/motif_curation/curated_30N_motif_pwmlist_20211014.rds" %>% readRDS()
#' seqs="/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned/B8_20_Mock__root__conc_0nM__time_2h__30n.gz"
#' motif_count <- motif_counts_(seqs,pwmlist)
motif_counts_ <- function(seqs,pwmlist){
  if((seqs %>% typeof())=="character" && length(seqs)==1){seqs=fjComm::getSeq_fqfachrFile(seqs)}
  dnabase <- seqs
  dnaseq <- DNAStringSet(dnabase)
  motif_ix_SummarizedExperiment <- motifmatchr::matchMotifs(pwmlist,dnabase,bg="even",out = "scores")
  motif_ct <- motifmatchr::motifCounts(motif_ix_SummarizedExperiment) %>% Matrix::colSums()
  motif_ffr <- motif_ct / length(dnaseq)

  #generate bk sequences with the same 1-nt frequency, match motif
  ntfreq <- Biostrings::alphabetFrequency(dnaseq, baseOnly=TRUE, collapse=TRUE) %>% prop.table()
  seqnons <-universalmotif::create_sequences(seqnum = length(dnaseq), seqlen = width(dnaseq) %>% tabulate() %>% which.max(),
                              freqs = ntfreq[c("A","C","G","T")])
  motif_ix_SummarizedExperiment_nomr <- motifmatchr::matchMotifs(pwmlist,seqnons,bg="even",out = "scores")
  motif_ct_nomr <- motifmatchr::motifCounts(motif_ix_SummarizedExperiment_nomr) %>% Matrix::colSums()
  motif_ffr_nomr <- motif_ct_nomr / length(seqnons)

  motif_deno <- motif_ffr - motif_ffr_nomr
  fil_motif <- motif_deno <= 0
  motif_deno[fil_motif] <- 0.000001 # adjust counts to 0 if less than 0
  motif_deno
}


#' plot motif enrichment of a target file or sequences
#'
#' @param seqs the sequence file, or a character of sequences
#' @param pwmlist a PWMatrixList containing the pwm of motifs
#'
#' @return a ggplot of the enrichment of each motif
#' @export
#'
#' @examples
#' pwmlist <- readRDS("/wrk/zhu/test/motif_curation/curated_30N_motif_pwmlist_20211014.rds")
#' seqs="/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned/B8_20_Mock__root__conc_0nM__time_2h__30n.gz"
#' tt= plot_motif_enrichment(seqs,pwmlist)
plot_motif_enrichment<-function(seqs,pwmlist,threads=1)
{
  # figure out the names of the motif
  parallelParam= BiocParallel::MulticoreParam(workers = threads)
  pwmlist <- BiocParallel::bplapply(pwmlist,function(x){motif_get_name(x)},BPPARAM = parallelParam) %>% suppressMessages()
          pwmlist_tmp <- TFBSTools::PWMatrixList()
            for(i in 1:length(pwmlist)){
              pwmlist_tmp[[i]] <- pwmlist[[i]]
            }
  pwmlist= pwmlist_tmp; rm("pwmlist_tmp")
  tfnames <- lapply(pwmlist,function(x)x@name) %>% unlist()
  motif_freqs <-motif_counts_(seqs,pwmlist) %>% suppressMessages()

  # # re-order with alphabetic
  #     order_=tfnames %>% order()
  #   tfnames=tfnames[order_]
  #   motif_freqs=motif_freqs[order_]
  #   pwmlist=pwmlist[order_]

    levs <- unique(tfnames)
    tfnames <- factor(tfnames,levels=levs)

  plot=ggplot2::ggplot() + geom_col(aes(tfnames,motif_freqs,fill=factor(tfnames)),width = 0.5,alpha=max(motif_freqs)*10/0.35)+
    coord_flip()+ scale_y_continuous(expand=c(0,0))+ fjComm::gg_theme_Publication() +
    ylab("Motif frequency")+xlab("")+theme(axis.ticks = element_blank())+ scale_fill_discrete(guide="none")


  motif_pics=BiocParallel::bplapply(1:length(pwmlist),
                                  function(x){
                                    t_pwmlist <- universalmotif::convert_motifs(pwmlist[[x]])
                                    motif_pic=suppressMessages(universalmotif::view_motifs(t_pwmlist,show.positions = F)+ylab("")+scale_y_continuous() +theme(axis.ticks.y =element_blank(),axis.line.y = element_blank(),axis.text.y = element_blank())  )
                                    ggplotGrob(motif_pic)
                                    },BPPARAM = BiocParallel::MulticoreParam(workers = threads))

  # construct plotting texts, then evaluate to generate the plot
            pwmlist_tmp= universalmotif::convert_motifs(pwmlist)
      img.name=  pwmlist_tmp %>% map(~.@name) %>% unlist()
        if(!identical(img.name,tfnames %>% as.character())) stop("TF ordering problem")
      width <- pwmlist_tmp %>% map(~.@consensus) %>% str_count() / 1.5 * 4
      # height <- pwmlist_tmp %>% map(~.@icscore) %>% unlist() %>% {. / 5 *4} %>% max(1.5)
    # labels_motif <- paste0("<img src='", motif_pics,  "' width=","'",width,"'" ," ","height='",height,"'"," /><br>*", img.name,"*")

      len=length(motif_pics)
      ymax_=max(motif_freqs);ymin_=ymax_-width/max(width)*ymax_/2
      anno_texts=glue::glue("annotation_custom(motif_pics[[{1:len}]],ymin = {ymin_},ymax = {ymax_},xmin = {1:len-1},xmax = {1:len+0.5})") %>% paste0(collapse = "+")
      plot_text=paste0("plot+ ",anno_texts)
  plot_w_motif=eval(parse(text = plot_text))

  plot_w_motif
}



# library(hbmcbioR)
# pwmlist <- readRDS("/wrk/zhu/test/motif_curation/curated_30N_motif_pwmlist_20211014.rds")
# seqs="/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned/B8_20_Mock__root__conc_0nM__time_2h__30n.gz"
# tt= plot_motif_enrichment(seqs,pwmlist,threads = 30)
# # fjComm::gg_save_pdf(tt,25,30,path = "/wrk/zhu/test/motif/",filename = "p_pic.pdf")

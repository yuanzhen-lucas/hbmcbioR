#' identify gave consensus or universalmotif motif to get TF html file
#'
#' @param pcm a motif of universalmotif format (PWM PCM or PPM)
#' @param n   the result number you want to get,defualt is 5
#' @param by  if consensus,use "by="consensus"".if universalmotif,use "by="pcm""
#'
#' @return
#' @export
#' @exaple
#'my_html <- motif_database_html("KDKDKDKAATCRRBB",by="consensus")
#'
#'motifdatabase=universalmotif::read_motifs("/wrk/yuanzhen/motifdb/motif_database")
#'pcm <- motifdatabase[[1]]
#'my_html <- motif_database_html(pcm,by="pcm")
motif_database_html <- function(pcm,n=5,by="pcm"){
  motifdatabase <- read_motifs("/wrk/yuanzhen/motifdb/motif_database")
  pwm_library <- read_rds("/wrk/yuanzhen/motifdb/pwm_library.Rds")
  if(by=="pcm"){
    pcm <- pcm
  }else if(by=="consensus"){
    pcm <- universalmotif::create_motif(pcm,name = "your_give_consensus",alphabet = "DNA")
  }else{
    stop("you parameter is not valid,please use pcm or consensus")
  }

  motiflogo <- universalmotif::convert_motifs(pcm, "TFBSTools-PWMatrix")
  # find the most similar motif to our motif
  pwm_sim <- TFBSTools::PWMSimilarity(pwm_library, motiflogo, method = 'Pearson')
  pwm_sim_rc <- TFBSTools::PWMSimilarity(pwm_library,convert_motifs( universalmotif::motif_rc(pcm),"TFBSTools-PWMatrix"), method = 'Pearson')
  pwm_sim <- round(pmax(pwm_sim,pwm_sim_rc),3)
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
  x_n <- head(pwm_library_dt,n)
  ge <- readr::read_rds("/wrk/yuanzhen/motifdb/motif_gene.Rds")
  name_ <- lapply(motifdatabase,function(x)x["name"]) %>% unlist()
  dtlogo <- .da_sim(x_n=x_n,pcm=pcm,pwm_library=pwm_library,motifdatabase=motifdatabase,name_=name_,ge=ge)
  dtlogoresult <- .addmotifLogorow(dtlogo,name_m=pcm["name"])
  DT::datatable(dtlogoresult,
                escape = FALSE, # To show the logo
                filter="top", options=list(pageLength=10))

}



#' find motif information in database
#'
#' @param x_n a data.frame include tf_name,similarity,and the location in database
#' @param pcm you give consensus or universalmotif format motif
#' @param pwm_library a PWMartrixList of database motif
#' @param motifdatabase a universalmotif format list of database motif
#' @param name_ the number of location
#' @param ge the text file of gene_id and TF_name,famliy or other informations
#'
#'
#' @return
#' @export
#'
#' @examples
#' @seealso [motif_database_html()]
.da_sim <- function(x_n,pcm,pwm_library,motifdatabase,name_,ge){
  dt_total <- data.frame()
  na_li <- x_n$col
  for(i in 1:length(x_n$col)){
    pwmnub <- na_li[[i]] %>% as.numeric()
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
      method <- name$method
      TF_name <- name$TF_name
      TF_family <- name$TF_family
      PMID <- name$PMID
    }else if(nrow(name)>1){
      method <- str_c(unique(name$method),collapse = "/")
      TF_name <- name$TF_name[[1]]
      TF_family <- name$TF_family[[1]]
      PMID <- str_c(unique(name$PMID),collapse = "/")
    }else{
      method <- "NO_source_method"
      TF_name <- "NO_source_name"
      TF_family <- "NO_source_family"
      PMID <- "NO_source_PMID"
    }

    motif_pic <- universalmotif::view_motifs(c(know_motif,pcm))
    ggplot2::ggsave(filename = paste0(TF_name,"_",x_n$similarity[[i]],".png"), plot = motif_pic,"/var/www/html/logos/logoraw/",device = "png")
    dtlogo <- data.frame(motifTF_name = TF_name,
                         motifTF_family = TF_family,
                         motifTF_similarity = x_n$similarity[[i]],
                         motifTF_TFlocus = tf_locus,
                         motifTF_method=method,
                         PMID=PMID,
                         database_source=database_source)
    dt_total <- rbind(dt_total,dtlogo)
  }
  dt_total
}




#' make html index
#'
#' @param motifEnrDT a dataframe include motifTF_name,motifTF_family,motifTF_similarity,and so on
#' @param addHTML add picture index ,default is TRUE
#' @param motifcolumname the column name of motif in motifEnrDT
#' @param motifsimi the column name of similarity in motifEnrDT
#' @param name_m the motif name
#'
#'
#'
#' @return
#' @export
#'
#' @examples
#' @seealso [motif_database_html()]
.addmotifLogorow <- function(motifEnrDT,addHTML=TRUE, motifcolumname = "motifTF_name",motifsimi="motifTF_similarity",name_m=pcm["name"])
{
  isNA <- which(motifEnrDT[[motifcolumname]]=="")
  logos <- paste("http://59.79.233.205/", "/","logos/logoraw/", motifEnrDT[[motifcolumname]],"_",motifEnrDT[[motifsimi]],".png", sep="")

  if(addHTML)
  {
    logos <- paste('<img src="', logos,
                   '" height="52" alt="',
                   motifEnrDT[[motifcolumname]], '"></img>', sep="")
  }
  logos[isNA] <- ""
  data.table::data.table(motifrowlogo=logos, motifEnrDT)
}

######## motif utilities ##########
#' Motif Discovery with Autoseed
#'
#'
#' @param sig_file  a file containing the reads for motif discovery, can be either a fastq, fasta, or plain text file
#' @param outdir  the path to store output data (seeds and svg visualizations)
#' @param bk_file the background file for motif discovery, leave blank to use random sequences
#' @param local_max_cutoff  the sensitivity of motif discovery (it tells at least how many hits of the seed consensus should exist in the sig_file, in order to generate a PFM output for this seed )
#'
#' @return
#' a list of pfms discovered from the sig_file, can use ggseqlogo for further visualization
#' @export
#'
#' @examples
#' sig_file="/var/www/html/selex/data/sequence/random40ND8.44.ATI4.seedling.root.20101n.1..30n4u_sig.seq"
#' result=autoseed_motif_disc(sig_file,"/wrk/zhu/test/autoseed_test",local_max_cutoff = 20)
autoseed_motif_disc<-function(sig_file,outdir=NULL,bk_file="/var/www/html/selex/data/sequence/Empty.seq",local_max_cutoff=40)
{
  wd=getwd()
  tmp_dir=tempfile(); dir.create(tmp_dir,recursive = TRUE)
    if(is.null(outdir)) outdir=tmp_dir #discard figure output if not needed
  tmp_file=tempfile(tmpdir = tmp_dir)
  setwd(tmp_dir)
        # adjust ligand length to 40bp for motif discovery
          tmp_seq_file=paste0(tmp_dir,"/",basename(sig_file)) %>% stringr::str_replace(".gz$",".txt")
          .length_to_40bp(sig_file,tmp_seq_file)
          sig_file=tmp_seq_file

  if(!dir.exists(outdir))dir.create(outdir,recursive = TRUE)
  cmd="/var/www/html/selex/totalautoseed4_1 -40N {bk_file} {sig_file} 1 8 10 0.35 - {local_max_cutoff} 100 {outdir} >{tmp_file}" %>% glue::glue()
  system(cmd)
  result=read_lines(tmp_file)
  unlink(tmp_dir,recursive = TRUE)
  setwd(wd)
  motif_data=result[result %>% str_detect("All Hits")]
    if(length(motif_data)==0) return(list(matrix(0.25,4,6)))
  .callbk=function(x,ind){x[3:6] %>% read_tsv(col_names = F) %>% .[,-c(1,2)] %>% as.matrix()}
  motif_list=read_lines_chunked(file = motif_data, callback= ListCallback$new(.callbk), chunk_size = 6)
  motif_list %>% map(.pfm_trim)
}









########  aux functions ########
#' adjust sequence length to 40bp in order to use autoseed for motif discovery
#'
#' @param in_file
#' @param out_file
#' @param segment_x3
#'
#' @return
#' @export
#'
#' @examples
.length_to_40bp<-function(in_file,out_file,segment_x3=TRUE)
{
  filename=in_file %>% basename()
  seqs= fjComm::getSeq_fqfachrFile(in_file)
  seqlens=seqs %>% nchar()

  tt= seqs %>% str_sub(1,40)
  Ns=40-nchar(tt)
  topaste=stringi::stri_rand_strings(length(Ns),Ns,"[ACGT]")
  tt=paste0(tt,topaste)

  if(segment_x3){
    # take end if ligand length >40
    # take middle if ligand length>80
    tt_e=seqs[seqlens>40] %>% str_sub(-40,seqlens[seqlens>40])
    m_start=ceiling(seqlens[seqlens>(40*2)]/2)-ceiling(40/2); m_end=m_start+40-1
    tt_m=seqs[seqlens>(40*2)] %>% str_sub(m_start,m_end)

    tt=c(tt,tt_e,tt_m)
  }

  if (tt %>% nchar %>% table() %>% length() %>% {.>1}) stop("ligands are not of the same lengths -->{filename}" %>% glue::glue())

  write_lines(tt,out_file)
  out_file
}


#' trim low IC bases flanking a pfm
#'
#' @param pfm
#' @param include_flank_cutoff
#'
#' @return
#' @export
#'
#' @examples
.pfm_trim<-function(pfm,include_flank_cutoff=2)
  # 1. include flanking positions where the ratio between the most and least frequent bases was > 2 (include_flank_cutoff)
{
  freq_pfm=pfm %>% sweep(2, 0.001, "+") %>% sweep(2,colSums(.),FUN = "/")
  keepfilter=(colMaxs(freq_pfm)/(colMins(freq_pfm)+0.01))>=include_flank_cutoff
  if(is.null(keepfilter) || keepfilter == '') return()
  if(all(keepfilter == FALSE)) return()
  keepfilter[which(keepfilter) %>% {min(.):max(.)}]=TRUE
  if(is.null(keepfilter) || keepfilter == '') return()
  pfm[,keepfilter]

}



####### test code ######

# fjComm::clear_()
# sig_file="/var/www/html/selex/data/sequence/random40ND8.44.ATI4.seedling.root.20101n.1..30n4u_sig.seq"
# result=autoseed_motif_disc(sig_file,"/wrk/zhu/test/autoseed_test",local_max_cutoff = 20)
# re=result %>% ggseqlogo_lab_list()

#
# autoseed_motif_disc<-function(sig_file,outdir=NULL,bk_file="/var/www/html/selex/data/sequence/Empty.seq",local_max_cutoff=40)
# {
#   # browser()
#   wd=getwd()
#   tmp_dir=tempfile(); dir.create(tmp_dir,recursive = TRUE)
#   if(is.null(outdir)) outdir=tmp_dir #discard figure output if not needed
#   tmp_file=tempfile(tmpdir = tmp_dir)
#   setwd(tmp_dir)
#   # adjust ligand length to 40bp for motif discovery
#   tmp_seq_file=paste0(tmp_dir,"/",basename(sig_file)) %>% stringr::str_replace(".gz$",".txt")
#   .length_to_40bp(sig_file,tmp_seq_file)
#   sig_file=tmp_seq_file
#
#   if(!dir.exists(outdir))dir.create(outdir,recursive = TRUE)
#   cmd="/var/www/html/selex/totalautoseed4_1 -40N {bk_file} {sig_file} 1 8 10 0.35 - {local_max_cutoff} 100 {outdir} >{tmp_file}" %>% glue::glue()
#   system(cmd)
#   result=read_lines(tmp_file)
#   unlink(tmp_dir,recursive = TRUE)
#   setwd(wd)
#   motif_data=result[result %>% str_detect("All Hits")]
#   .callbk=function(x,ind){x[3:6] %>% read_tsv(col_names = F) %>% .[,-c(1,2)] %>% as.matrix()}
#   motif_list=read_lines_chunked(file = motif_data, callback= ListCallback$new(.callbk), chunk_size = 6)
#
#   for(i in 1:length(motif_list))
#   {
#     print(paste0("cycle",i))
#     # if(i==5) browser()
#     .pfm_trim(motif_list[[i]])
#     }
#   motif_list %>% map(.pfm_trim)
# }
#
# .pfm_trim<-function(pfm,include_flank_cutoff=2)
#   # 1. include flanking positions where the ratio between the most and least frequent bases was > 2 (include_flank_cutoff)
# {
#   # browser()
#   freq_pfm=pfm %>% sweep(2,colSums(.),FUN = "/")
#   print(pfm)
#   keepfilter=(colMaxs(freq_pfm)/(colMins(freq_pfm)+0.01))>=include_flank_cutoff
#   keepfilter[which(keepfilter) %>% {min(.):max(.)}]=TRUE
#   pfm[,keepfilter]
# }
#

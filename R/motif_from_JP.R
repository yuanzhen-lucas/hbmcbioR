#' make a pcm to the most similar known TFs motif name
#' @param pcm  include your motif pcm
#' @param species your interested species(the JASPAR included),default is Arabidopsis thaliana
#' @param n the number of similar motif name,default is 1
#'
#'
#' @example
#' my_know <- motif_from_JP(root_motifs[[1]])


motif_from_JP <- function(pcm,species='Arabidopsis thaliana',n=1){
#to get all of the TFs PWM from JASPAR2020
  pwm_library <- getMatrixSet(
    JASPAR2020,
    opts=list(
      species =  species,
      matrixtype = 'PWM'
    ))

    motiflogo <- convert_motifs(pcm, "TFBSTools-PWMatrix")
    # find the most similar motif to our motif
    pwm_sim <- PWMSimilarity(pwm_library, motiflogo, method = 'Pearson')
    pwm_sim_rc <- PWMSimilarity(pwm_library, universalmotif::motif_rc(motiflogo), method = 'Pearson')
    pwm_sim <- pmax(pwm_sim,pwm_sim_rc)
    # extract the motif names from the pwm library
    pwm_library_list = lapply(pwm_library, function(x){
      data.frame(ID = ID(x))})

    # combine the list into one data frame
    pwm_library_dt = dplyr::bind_rows(pwm_library_list)

    # fetch the similarity of each motif to our unknown motif
    pwm_library_dt$similarity = pwm_sim[pwm_library_dt$ID]

    # find the most similar motif in the library
    pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]
    x_n <- head(pwm_library_dt,n)
    pfm_all <- sapply(x_n$ID,function(x){TFBSTools::getMatrixByID(JASPAR2020, ID = x)})
    lapply(pfm_all, function(x)x@name)

}


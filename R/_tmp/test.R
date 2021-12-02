
filedress <- list.files("/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned")
filename <- paste0("/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned", "/", filedress)

for(i in 1:length(filedress)){
  if (filedress[i]!="B11_23_Mock101bp__root__conc_0nM__time_30min__101n.gz")next()
  varnames <- gsub('?.gz','',filedress[i])
  cat("Now will deal with the sample", "", filedress[i], "\n")
  dt <- autoseed_motif_disc(filename[i])
  cat("the file is", "",filename[i],"\n")
  write_rds(dt, paste0("/wrk/yuanzhen/motifdicr/",varnames,".Rds"))
  cat("we finished the sample", "",filedress[i], "\n")

}


library(BiocParallel)
filedress <- list.files("/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned")
filename <- paste0("/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned", "/", filedress)

f <- function(x){
  autoseed_motif_disc(x)
}
bplist <- bplapply(filename, f, BPPARAM = MulticoreParam(workers=30,progressbar=TRUE))
dt <- bplist(unlist)

my_packages <- c("universalmotif", "cowplot", "ggtree", "dplyr", "ggplot2", "TFBSTools", "JASPAR2020", "GenomicRanges","SummarizedExperiment", "Biostrings", "stats", "plotly", "htmlwidgets", "motifmatchr", "ShortRead","BiocParallel","DT","stringr","tidyverse","purrr")
lapply(my_packages, require, character.only = TRUE)

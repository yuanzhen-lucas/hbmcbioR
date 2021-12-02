outpath="analysis/2_0_enrich_plot/"; dir.create(outpath,recursive = TRUE,showWarnings = F)
bamfiles= sample_info$bam_files
pacman::p_load(GenomicRanges,GenomicAlignments,GenomicFeatures,Rsamtools,EnrichedHeatmap,GenomicFeatures)
txdb=makeTxDbFromGFF(gffFile)
# if (is.null(gffFile)) {
#   pacman::p_load(TxDb.Athaliana.BioMart.plantsmart51)
#   txdb=TxDb.Athaliana.BioMart.plantsmart51
# }

vis_width=vis_range*2+1
tras <- transcripts(txdb)
tss=promoters(tras,upstream = 0,downstream = 1)
tss=resize(tss,vis_width,fix = "center") %>% trim()
tss=tss[width(tss)==vis_width]
rngs=ranges(tss)
start_=IRanges::start(rngs)
filter_=start_>0
tss=tss[filter_]

## remove 属于线粒体与叶绿体的区间
seqlevels(tss) %<>% str_replace("Chr","chr")
seqlevels(tss,pruning.mode="coarse")=seqlevels(tss)[seqlevels(tss) %>% str_detect("[0-9]$")]


# lapply(seq_along(bamfiles), function(i){
parallel(seq_along(bamfiles), function(i){
  print(glue::glue("step{step} processing file: ") %>% paste0( bamfiles[i] %>% basename()))
  cov=coverage(bamfiles[i])
  for (j in 1:length(cov))
  {
    cov[[j]]=c(cov[[j]],Rle(0L,vis_width*2))
  }
  # ajust names according to names in tss
  names(cov) %<>% str_replace("Chr","chr")
  cov=cov[tss]
  cov %<>% unlist(use.names = F) %>% as.integer()
  cov=matrix(cov,ncol=vis_width,byrow = T)
  # flip genes on the reverse strand
  tss_minus=(strand(tss)=="-") %>% as.logical()
  cov[tss_minus,]= cov[tss_minus,ncol(cov):1]
  cov_acc=colSums(cov)
  cov_plot=ggplot()+geom_line(aes(-vis_range:vis_range,cov_acc))+scale_x_continuous(expand = c(0,0))+xlab("TSS")+ylab("MNase-seq coverage")
  gg_save_png(cov_plot,6,5,filename = outpath %>% paste0("plot/",fileNames[i],"_TSS"))
  gg_save_pdf(cov_plot,6,5,filename = outpath %>% paste0("plot/",fileNames[i],"_TSS"))
},workers = threads)
# })














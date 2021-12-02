uniquemotif <- function(indir,outdir,webdir="/var/www/html/logos",ic_threshold=1,icscore_threshold = 3,iclength_treshold=0.3,threshold=0.95,thread = 30){

  cat("we will deal with your data by four steps","\n")
  cat("1. de novo motif discovry","\n")
  cat("2. find the unique motif","\n")
  cat("3. detect the TF name of the unique motif from JASPAR2020","\n")
  cat("4. visualize the results","\n")


  ###Now we will the first steps to de novo motif discovry
  cat("Now we will the first steps to de novo motif discovry","\n")

  filedressdisc <- list.files(indir, pattern = ".gz")
  filenamedisc <- paste0(indir, "/", filedressdisc)

  motifdicr <- paste0(outdir,"/","motifdicr")
  dir.create(motifdicr)


  f <- function(x){
    autoseed_motif_disc(x)
  }

  dtdicr <- bplapply(filenamedisc, f, BPPARAM = MulticoreParam(workers=30,progressbar=TRUE))

  write_rds(dtdicr,paste0(outdir,"/","total_matrix",".Rds"))



  for(c in 1:length(dtdicr)){
    for(s in 1:length(dtdicr[[c]])){
      write_rds(dtdicr[[c]][[s]], paste0(motifdicr,"/",c,"_",s,".Rds"))
    }
  }

  motifrawdir <- paste0(outdir,"/","motifrawdir")
  dir.create(motifrawdir)

  sum_comse <- function(con){
    nub <- con %>% str_split("") %>%.[[1]] %>% as.data.frame()
    colnames(nub) <- "base"
    nub$nuber <- 1:length(nub$base)
    for(i in 1:length(nub$base)){
      if(nub$base[i] == "A"){
        nub[i,3] <- 1
      }else if(nub$base[i] == "C"){
        nub[i,3] <- 2
      }else if(nub$base[i] == "G"){
        nub[i,3] <- 3
      }else if(nub$base[i] == "T"){
        nub[i,3] <- 4
      }else if(nub$base[i] == "N"){
        nub[i,3] <- 5
      }else if(nub$base[i] == "Y"){
        nub[i,3] <- 6
      }else if(nub$base[i] == "W"){
        nub[i,3] <- 7
      }else if(nub$base[i] == "S"){
        nub[i,3] <- 8
      }else if(nub$base[i] == "M"){
        nub[i,3] <- 9
      }else if(nub$base[i] == "K"){
        nub[i,3] <- 10
      }else if(nub$base[i] == "H"){
        nub[i,3] <- 11
      }else if(nub$base[i] == "B"){
        nub[i,3] <- 12
      }else if(nub$base[i] == "V"){
        nub[i,3] <- 13
      }else if(nub$base[i] == "D"){
        nub[i,3] <- 14
      }else if(nub$base[i] == "R"){
        nub[i,3] <- 15
      }else{
        next()
      }
    }
    colnames(nub) <- c("base","nuber","value")
    nub$multi <- nub$nuber * nub$value
    sum(nub$multi)
  }


  matrixlist <- list()
  filedressmatrix <- list.files(motifdicr)
  filenamematrix <- paste0(motifdicr, "/", filedressmatrix)

  for(v in 1:length(filedressmatrix)){
    n_name <- gsub('?.Rds','',filedressmatrix[v])
    a <- read_rds(filenamematrix[v])
    if(is.null(a)) next()
    if(is.null(nrow(a))) next()
    b <- as.matrix(a)
    colnames(b) <- NULL
    rownames(b) <- c("A", "C", "G", "T")
    motifyz <- convert_motifs(b)
    motifyzrc <- motif_rc(motifyz)
    nub <- sum_comse(motifyz["consensus"])
    nubrc <- sum_comse(motifyzrc["consensus"])
    if(nub > nubrc){
      motifyzrc["name"] <- n_name
      write_motifs(motifyzrc, paste0(motifrawdir, "/", "motif", n_name))
    }else{
      motifyz["name"] <- n_name
      write_motifs(motifyz, paste0(motifrawdir, "/", "motif", n_name))
    }
  }

  cat("2. find the unique motif","\n")
  IC_calc<-function(consensus,pseudo_freq=0.0000000000001)
  {
    letters=consensus %>% str_split("") %>% .[[1]]
    letter_freqs= letters %>% table() %>% prop.table()
    for (letter in c("A","C","G","T"))
    {
      if (is.na(letter_freqs[letter])) letter_freqs[letter]=pseudo_freq
    }
    ic=2-sum(-letter_freqs*(log2(letter_freqs)))
    ic
  }

  filedress <- list.files(motifrawdir)
  filename <- paste0(motifrawdir, "/",filedress)

  for(y in 1:length(filename)){
    motiftrim <- read_motifs(filename[y])
    if(str_detect(motiftrim["consensus"],"^N+[^ATCGRYWSMKHBVD]$")){
      cat("we remove motif", filename[y],"\n","because of repeat N consensus:","",motiftrim["consensus"],"\n")
      unlink(filename[y])
    }else if(str_detect(motiftrim["consensus"],"^N*[[RYWSMKHBVD]+N*[RYWSMKHBVD]+]+N*$")){
      cat("we remove motif", filename[y],"\n","because of unreliable consensus:","",motiftrim["consensus"],"\n")
      unlink(filename[y])
    }else if(IC_calc(motiftrim["consensus"]) > ic_threshold){
      cat("we remove motif", filename[y],"\n","because of IC too hight","",motiftrim["consensus"],"\n")
    }else{
      y=y+1
    }
  }


  for (n in 1:length(filename)){

    if(file.exists(filename[n])){
      t1 <- read_motifs(filename[n])
      t1rc <- motif_rc(t1)
      pwm_1 <- convert_motifs(t1, "TFBSTools-PWMatrix")
      pwm_1rc <- convert_motifs(t1rc, "TFBSTools-PWMatrix")
      cat("we deal with motif",filename[n],"just step one","\n")
      p <- n+1
      leqm <- length(filename)
      for (m in p:leqm){
        if(file.exists(filename[m])){
          t2 <- read_motifs(filename[m])
          pwm_2 <- convert_motifs(t2, "TFBSTools-PWMatrix")
          if(((PWMSimilarity(pwm_1, pwm_2, method="Pearson")) > threshold)){
            if(str_count(t1["consensus"],"N") > str_count(t2["consensus"],"N")){
              cat("we remove motif", filename[m],"it's consensus",t1["consensus"],"\n")
              unlink(filename[n])
            }else if(str_count(t1["consensus"],"N") < str_count(t2["consensus"],"N")){
              cat("we remove motif", filename[m],"it's consensus",t2["consensus"],"\n")
              unlink(filename[m])
            }else{
              cat("we remove motif", filename[m],"it's consensus",t2["consensus"],"\n")
              unlink(filename[m])
            }

          } else if((PWMSimilarity(pwm_1rc, pwm_2, method="Pearson")) > threshold) {

            if(str_count(t1rc["consensus"],"N") > str_count(t2["consensus"],"N")){
              cat("we remove motif", filename[m],"it's consensus",t1["consensus"],"\n")
              unlink(filename[n])
            }else if(str_count(t1rc["consensus"],"N") < str_count(t2["consensus"],"N")){
              cat("we remove motif", filename[m],"it's consensus",t2["consensus"],"\n")
              unlink(filename[m])
            }else{
              cat("we remove motif", filename[m],"it's consensus",t2["consensus"],"\n")
              unlink(filename[m])
            }
          }else{
            m=m+1
          }
        }
      }
    }else{
      n=n+1
    }
  }

  motifunique <- list()
  fileunique <- list.files(motifrawdir)
  filenameunique <- paste0(motifrawdir,"/",fileunique)
  for(h in 1:length(fileunique)){
    motifunique <- c(motifunique,list(read_motifs(filenameunique[h])))
  }
  motif_dt_unique <- to_df(motifunique)

  motif_dt_unique_1 <- mutate(motif_dt_unique,iclength = icscore / str_length(consensus))


  motiffilterdir <- paste0(outdir,"/","motiffilterdir")
  dir.create(motiffilterdir)

    motiffilter_dt <- filter(motif_dt_unique_1,icscore > icscore_threshold)
    motiffilter_dt_1 <- filter(motiffilter_dt,iclength > iclength_treshold)
    name_filter <- motiffilter_dt_1$name

  sapply(paste0("motif",name_filter),function(x){file.copy(paste(motifrawdir,x,sep="/"),motiffilterdir)})
  #mapply(file.copy, from=paste0(motifrawdir,"/","motif",name_filter), to=motiffilterdir)

  matrixlistlast <- list()
  filedressmatrixlast <- list.files(motiffilterdir)
  filenamematrixlast <- paste0(motiffilterdir, "/", filedressmatrixlast)
  basetest <- str_extract(basename(filenamematrixlast),"[0-9]+_[0-9]+")
  cholist  <- lapply(basetest, function(x){
    read_rds(paste0(motifdicr,"/",x,".Rds"))
  })

  fjComm::ggseqlogo_lab_list(cholist)


  ggseqlogo_save_A4_pdf <- function(mat_list, filePrefix="motifplots_p", motifs_per_page=60,...)
  {
    #name147=lapply(mat_list,function(x)attr(x,"name")) %>% unlist(); names(mat_list)=name147 %>% str_replace("_.*$","") %>% str_replace("ooo.*$","")
    mat_list %<>% purrr::map(~{set_rownames(as.matrix(.),c("A","C","G","T"))})
    length=length(mat_list); pages=(length/motifs_per_page) %>% as.integer()+1
    for (page in 1:pages)
    {
      start_=(motifs_per_page*(page-1)+1)
      end_=motifs_per_page*page; if(end_>length) end_=length
      # if(page==4)browser()
      rows=(((start_:end_) %>% length )/5) %>% base::ceiling()
      p147=mat_list[start_:end_] %>% ggseqlogo( ncol=5,...)+ scale_y_continuous(breaks = c(0,2)) +gg_theme_Publication()+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank() ) #theme(axis.text.x = element_blank(),title = )
      gg_save_pdf(p147,21,29.7/12*rows,filename = filePrefix %>% paste0(page))
    }
  }






  cat("we will detect the TF name of the unique motif from JASPAR2020","\n")

  pwm_library <- getMatrixSet(
    JASPAR2020,
    opts=list(
      species    = 'Arabidopsis thaliana',
      matrixtype = 'PWM'
    ))

  dtlogo <- data.frame()
  filedresslogo <- list.files(motiffilterdir)
  filenamelogo <- paste0(motiffilterdir,"/",filedresslogo)
  logoraw <- paste0(webdir,"/","logoraw")
  logoTF <- paste0(webdir,"/","logoTF")
  basewebname <- basename(webdir)
  dir.create(logoraw)
  dir.create(logoTF)


  for(d in 1:length(filedresslogo)){

    motiflogoraw <- read_motifs(filenamelogo[d])
    motiflogo <- convert_motifs(motiflogoraw, "TFBSTools-PWMatrix")
    # find the most similar motif to our motif
    pwm_sim <- PWMSimilarity(pwm_library, motiflogo, method = 'Pearson')

    # extract the motif names from the pwm library
    pwm_library_list = lapply(pwm_library, function(x){
      data.frame(ID = ID(x))})

    # combine the list into one data frame
    pwm_library_dt = dplyr::bind_rows(pwm_library_list)

    # fetch the similarity of each motif to our unknown motif
    pwm_library_dt$similarity = pwm_sim[pwm_library_dt$ID]

    # find the most similar motif in the library
    pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]

    name1 <- motiflogoraw["name"]
    name2 <- pwm_library_dt$ID[1]
    name3 <- pwm_library_dt$ID[2]
    name4 <- pwm_library_dt$ID[3]


    pfm1 <- TFBSTools::getMatrixByID(JASPAR2020, ID = name2)
    pfm2 <- TFBSTools::getMatrixByID(JASPAR2020, ID = name3)
    pfm3 <- TFBSTools::getMatrixByID(JASPAR2020, ID = name4)

    dtlogo <- rbind(dtlogo,data.frame(motifrowname = name1,
                                      motifTFname1 = pfm1@name,
                                      motifTFname2 = pfm2@name,
                                      motifTFname3 = pfm3@name))


    pngmotif <- view_motifs(motiflogoraw,use.type = "PPM", sort.positions = T)
    ggsave(filename = paste0(name1,".png"), plot = pngmotif,logoraw,device = "png")

    pngmotif1 <- view_motifs(pfm1,use.type = "PPM", sort.positions = T)
    ggsave(filename = paste0(name1,"_",pfm1@name,".png"), plot = pngmotif1,logoTF,device = "png")

    pngmotif2 <- view_motifs(pfm2,use.type = "PPM", sort.positions = T)
    ggsave(filename = paste0(name1,"_",pfm2@name,".png"), plot = pngmotif2,logoTF,device = "png")

    pngmotif3 <- view_motifs(pfm3,use.type = "PPM", sort.positions = T)
    ggsave(filename = paste0(name1,"_",pfm3@name,".png"), plot = pngmotif3,logoTF,device = "png")


  }
  write_rds(dtlogo,paste0(outdir,"/","dtlogo",".Rds"))

  cat("we will visualize the results","\n")

  addmotifLogotf1 <- function(motifEnrDT,addHTML=TRUE, motifcolumname = "motifrowname", motifTFname)
  {
    isNA <- which(motifEnrDT[[motifcolumname]]=="")

    logos <- paste("http://59.79.233.205/","/",basewebname, "/","logoTF/", motifEnrDT[[motifcolumname]],"_",motifEnrDT[[motifTFname]],".png", sep="")

    if(addHTML)
    {
      logos <- paste('<img src="', logos,
                     '" height="52" alt="',
                     motifEnrDT[[motifcolumname]], "_",motifEnrDT[[motifTFname]],'"></img>', sep="")
    }

    logos[isNA] <- ""

    if(motifTFname == "motifTFname3"){
      data.table::data.table(motifTFname3=logos, motifEnrDT)
    }else if(motifTFname == "motifTFname2"){
      data.table::data.table(motifTFname2=logos, motifEnrDT)
    }else{
      data.table::data.table(motifTFname1=logos, motifEnrDT)
    }

  }

  dtlogomotiftf3 <- addmotifLogotf1(dtlogo,motifTFname = "motifTFname3")
  dtlogomotiftf2 <- addmotifLogotf1(dtlogomotiftf3,motifTFname = "motifTFname2")
  dtlogomotiftf1 <- addmotifLogotf1(dtlogomotiftf2,motifTFname = "motifTFname1")

  addmotifLogorow <- function(motifEnrDT,addHTML=TRUE, motifcolumname = "motifrowname")
  {
    isNA <- which(motifEnrDT[[motifcolumname]]=="")
    logos <- paste("http://59.79.233.205/","/",basewebname, "/","logoraw/", motifEnrDT[[motifcolumname]],".png", sep="")

    if(addHTML)
    {
      logos <- paste('<img src="', logos,
                     '" height="52" alt="',
                     motifEnrDT[[motifcolumname]], '"></img>', sep="")
    }
    logos[isNA] <- ""

    data.table::data.table(motifrowlogo=logos, motifEnrDT)
  }

  dtlogoresult <- addmotifLogorow(dtlogomotiftf1)

  yz <- DT::datatable(dtlogoresult,
                      escape = FALSE, # To show the logo
                      filter="top", options=list(pageLength=10))

  htmlfile <- paste0(outdir,"/","result",".","html")

  DT::saveWidget(yz, htmlfile)


}


#' you can compare two fasta file include the enrichment of your interested motifs
#'
#'
#' @param fastaFile1  a file containing the reads must be fasta file
#' @param fastaFile2  a file containing the reads must be fasta file
#' @param motifdir the directory containing the motifs(universalmotif format)
#'
#'
#' @return
#'a histogram png file
#' @export
#'
#' @examples
#'  motif_list_sample("/wrk/yuanzhen/b6test/30.fa","/wrk/yuanzhen/b6test/101.fa","/wrk/yuanzhen/b6test/motiffilterdir")
#' my_png=motif_list_sample("/wrk/yuanzhen/b6test/30.fa","/wrk/yuanzhen/b6test/101.fa","/wrk/yuanzhen/b6test/motiffilterdir")
motif_list_sample <- function(fastaFile1,fastaFile2,motifdir){

  fileplotlist <- list.files(motifdir)
  fileplotlistname <- paste0(motifdir,"/",fileplotlist)
  pwmlist <- PWMatrixList()

  for(f in 1:length(fileplotlist)){
    pwmlist[[f]] <- convert_motifs(read_motifs(fileplotlistname[f]), class = "TFBSTools-PWMatrix")
  }

  # extract the motif names from the pwm library
  pwm_library_list = lapply(pwmlist, function(x){
    data.frame(name = x@name)
  })

  # combine the list into one data frame
  pwm_library_dt = dplyr::bind_rows(pwm_library_list)

  seqs1 <- readDNAStringSet(fastaFile1)
  motif_ix_SummarizedExperiment1 <- matchMotifs(pwmlist,seqs1)

  countmotif1 <- motifMatches(motif_ix_SummarizedExperiment1)
  dt11 <- Matrix::summary(countmotif1)
  countnub1 <- table(dt11$j)
  dt21 <- as.data.frame(countnub1)

  tfcount1 <- cbind(dt21, pwm_library_dt)

  dt_total1 <- data.frame(tfname = tfcount1$name,
                          normseq = tfcount1$Freq)
  dt_total1$seqlength <- length(seqs1)

  dt_total1$ratio <- apply(dt_total1,1,function(x){x[2] %>% as.integer / x[3] %>% as.integer})



  seqs <- readDNAStringSet(fastaFile2)
  motif_ix_SummarizedExperiment <- matchMotifs(pwmlist,seqs)

  countmotif <- motifMatches(motif_ix_SummarizedExperiment)
  dt1 <- Matrix::summary(countmotif)
  countnub <- table(dt1$j)
  dt2 <- as.data.frame(countnub)

  tfcount <- cbind(dt2, pwm_library_dt)

  tfcount$seqlength <- length(seqs)


  dt_total <- data.frame(tfname = tfcount$name,
                         normseq = tfcount$Freq)
  dt_total$seqlength <- length(seqs)
  dt_total$ratio <- apply(dt_total,1,function(x){x[2] %>% as.integer / x[3] %>% as.integer / 3})


  dt_total1$filename <- basename(fastaFile1)
  dt_total$filename <- basename(fastaFile2)
  dt_total1$num <- 2
  dt_total$num <- 1
  tdt <- rbind(dt_total,dt_total1)
  tdt <- arrange(tdt,tfname)
  yz <- ggplot(data = tdt, mapping = aes(x = num, y = ratio,fill = filename)) +
    geom_bar(stat = 'identity') + facet_grid(tfname~.)+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    fjComm::gg_theme_Publication() +
    coord_flip()
  yz

}




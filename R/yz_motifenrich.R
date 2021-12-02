#'split the total matrix file
#'
#'@param inputfile the totall motif matrix txt file
#'@param output the split results
#'
#'@export
#'@example
#'splitFile("totall.txt","outdir")


splitFile <- function(inputfile, output){
  # Create vector of packages
  my_packages <- c("universalmotif", "cowplot", "ggtree", "dplyr",
                   "ggplot2", "TFBSTools", "JASPAR2020")
  # Load multiple packages
  lapply(my_packages, require, character.only = TRUE)

  splitDir <- paste0(output, "/", "splitmotifdir")
  dir.create(splitDir)
  totalfile <- file(inputfile, open = "rt")
  rownb <- 4
  varnames <- gsub('?.txt','',inputfile)
  var_name <- paste(varnames, 1:1000, sep = "_")
  i <- 1
  while(TRUE){
    assign(var_name[i],value = readLines(totalfile,n = rownb,warn = FALSE))

    write.table(get(var_name[i]),
                paste0(splitDir, "/", var_name[i], ".txt"), row.names = F,col.names = F, quote = FALSE)
    if (length(get(var_name[i])) < rownb){
      unlink(paste0(splitDir, "/", var_name[i], ".txt"))
      break
    }
    else i <- i + 1
    }

}


#'chose the unique motif matrix
#'
#'@param splitmotifdir the dir include the split motif matrix file
#'@param outdir the dir remove the repeat motif matrix
#'
#'@export
#'@example
#'copMotif("outdir/splitmotifdir","outdir")
copMotif <- function(inputdir, outdir){
  # Create vector of packages
  #my_packages <- c("universalmotif", "cowplot", "ggtree", "dplyr",
  #"ggplot2", "TFBSTools", "JASPAR2020")
  # Load multiple packages
  #lapply(my_packages, require, character.only = TRUE)

  outdir <- "outdir"
  inputdir <- "outdir/splitmotifdir"
  ## make the PFM to universalmotif PFM
  uniqdir <- paste0(outdir, "/", "uniqdir")
  dir.create(uniqdir)
  txttempdir=paste0(outdir, "/", "txttempdir")
  dir.create(txttempdir)
  filedress <- list.files(inputdir)
  sapply(filedress,function(x){file.copy(paste(inputdir,x,sep="/"),txttempdir)})
  filedress <- list.files(txttempdir)
  filename <- paste0(txttempdir, "/", filedress)
  leq <- length(filedress)
  for (i in 1:leq){
    if(file.exists(filename[i]) == TRUE){
      n_name <- gsub('?.txt','',filedress[i])
      t1 <- read.delim(filename[i], header = FALSE)
      cat("we deal with",filename[i],"just step one","\n")
      n <- i+1
      leqm <- leq
      for (m in n:leqm){
        if(file.exists(filename[m]) == TRUE){
          t2 <- read.delim(filename[m], header = FALSE)
          if(setequal(t1, t2)== TRUE){
            unlink(filename[m])
            cat("we remove",filename[m],"\n")
          }else{
            m=m+1
          }
        }
      }
    }else{
      i=i+1
    }
  }
  filenewdress <- list.files(txttempdir)
  sapply(filenewdress,function(x){file.copy(paste(txttempdir,x,sep="/"),uniqdir)})
  unlink(txttempdir)
}



#'find the de novo motif with known TFs motif
#'
#'@param uniquedir the dir include the unique motif matrix file
#'@param outdir the motif logos
#'
#'@export
#'@example
#'compMotif("outdir/uniquemotif","outdir")

compMotif <- function(uniquedir, outdir){
  ## make the PFM to universalmotif PFM
  filedress <- list.files(uniquedir)
  motifdir=paste0(outdir, "/", "motifdir")
  dir.create(motifdir)

  pwm_library <- getMatrixSet(
    JASPAR2020,
    opts=list(
      species    = 'Arabidopsis thaliana',
      matrixtype = 'PWM'
    ))
  dt <- data.frame()

  for (i in 1:length(filedress)){
    filename <- paste0(uniquedir, "/", filedress)
    n_name <- gsub('?.txt','',filedress[i])
    a <- read.delim(filename[i], header=FALSE)
    b <- as.matrix(a)
    colnames(b) <- NULL
    rownames(b) <- c("A", "C", "G", "T")
    motif <- convert_motifs(b)
    motif["name"] <- n_name

    your_pwm <- convert_motifs(motif, "TFBSTools-PWMatrix")

    # find the most similar motif to our motif
    pwm_sim <- PWMSimilarity(pwm_library, your_pwm, method = 'Pearson')

    # extract the motif names from the pwm library
    pwm_library_list = lapply(pwm_library, function(x){
      data.frame(ID = ID(x), name = name(x))
    })

    # combine the list into one data frame
    pwm_library_dt = dplyr::bind_rows(pwm_library_list)

    # fetch the similarity of each motif to our unknown motif
    pwm_library_dt$similarity = pwm_sim[pwm_library_dt$ID]

    # find the most similar motif in the library
    pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]

    resultch <- head(pwm_library_dt, 1)

    name1 <- motif["name"]
    name2 <- pwm_library_dt$name[1]
    simi <- pwm_library_dt$similarity[1]

    motif["name"] <- paste0(name1, "_", name2)

    dt <- rbind(dt, data.frame(motifrowname = name1,
                               motifTFname = name2,
                               similarity = simi))


    pngmotif <- view_motifs(motif,use.type = "PPM", sort.positions = T)
    pngmotifTF <- view_motifs(pfmTF,use.type = "PPM", sort.positions = T)

    ggsave(filename = paste0(n_name,".png"), plot = pngmotif,
           path = "outdir/motifrow")
    ggsave(filename = paste0(n_name,".png"), plot = pngmotifTF,
           path = "outdir/motifTF")
    write_motifs(motif, paste0(motifdir, "/", "motif", n_name))

  }


  write.table(dt,paste0(outdir, "/", "motif", ".", "txt"), quote = F)
}


#'drow a fold changes and -log10 picture
#'
#'@param sequencefile your raw data sequence file must fasta file
#'@param prefix name of the output
#'@param backgroundseq background sequence also fasta file,and if NULL the file will  automatic generation
#'@export




plotEnrich <- function(sequencefile, prefix, backgroundseq=NULL){

  pwm_library <- getMatrixSet(JASPAR2020,
                              opts=list(species  = 'Arabidopsis thaliana',
                                        matrixtype = 'PWM'))

  start <- (proc.time())[3][[1]]

  seqs <- readDNAStringSet(sequencefile)

  motif_ix_SummarizedExperiment <- matchMotifs(pwm_library,seqs)

  countmotif <- motifMatches(motif_ix_SummarizedExperiment)
  dt1 <- Matrix::summary(countmotif)
  countnub <- table(dt1$j)
  dt2 <- as.data.frame(countnub)
  # extract the motif names from the pwm library
  pwm_library_list = lapply(pwm_library, function(x){
    data.frame(ID = ID(x), name = name(x))
  })

  # combine the list into one data frame
  pwm_library_dt = dplyr::bind_rows(pwm_library_list)

  tfcount <- cbind(dt2, pwm_library_dt)

  tfcount$seqlength <- length(seqs)



  if (is.null(backgroundseq) == T){
    seqnons <- create_sequences(seqnum = 5201314, seqlen = 30,
                                freqs = c(A=0.25, C=0.25, G=0.25, T=0.25))
    cat("we will produce random DNA sequence as background sequence,",
        "\n",
        "sequences number is 5201314,sequences length is 30,",
        "\n")
  } else{
    seqnons <- readDNAStringSet(backgroundseq)
    cat("we will deal with the background sequences as", backgroundseq, "\n")
  }

  motifnon <- matchMotifs(pwm_library,seqnons)
  countmotifnon <- motifMatches(motifnon)
  dt1non <- Matrix::summary(countmotifnon)
  countnubnon <- table(dt1non$j)
  dt2non <- as.data.frame(countnubnon)
  tfcountnon <- cbind(dt2non, pwm_library_dt)
  tfcountnon$seqnonlength <- length(seqnons)
  tfcountnon$seqslength <- length(seqs)
  normseqnon <- apply(tfcountnon,1,function(x){round(x[2] %>% as.numeric * x[6] %>% as.numeric / x[5] %>% as.numeric)})
  tfcountnon$normseqnon <- normseqnon

  dt_total <- data.frame(tfname = tfcount$name,
                         normseq = tfcount$Freq,
                         normseqnon = tfcountnon$normseqnon)
  pv <- apply(dt_total,1,function(x){poisson.test(x[[2]] %>% as.integer, x[[3]] %>% as.integer)})
  pvunlist <- lapply(pv,function(x){x$p.value}) %>% unlist()
  logp <- -log10(pvunlist)
  dt_total$logp <- logp
  dt_total$logp[which(logp > 300)] <- 300

  flod <- apply(dt_total,1,function(x){x[2] %>% as.integer / x[3] %>% as.integer})
  dt_total$flod <- flod
  dt_total$flod[which(flod > 10)] <- 10

  ggpic <- ggplot(data=dt_total, aes(x=logp, y=flod,color = tfname)) +
    geom_point(size=1) +
    ylim(0, 10) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(title="Known motif enrichment analysis of ATI data from sample l1",
         x="-log(p)", y="Fold change")

  ggplotmotif <- ggplotly(ggpic)
  filenames <- paste0(prefix,".html")
  saveWidget(ggplotmotif, filenames, selfcontained = F, libdir = "lib")
  end <- (proc.time())[3][[1]]
  cat("Aha, I quickly completed your task", "\n", "time = ", (end-start), "s", "\n")
}



#' make a html file
#'
#' @param motifEnrDT is a dataframe include motifTFname,motifrowname,similarity fuction:compMotif
#'
#' @param motifcolumname you can chose which to add the logo ,you must build a index
#' @export
#' @example
#' addmotifLogo(motifEnrDT,motifcolumname = "motifrowname")
addmotifLogo <- function(motifEnrDT,addHTML=TRUE, motifcolumname = "motifrowname")
{
  isNA <- which(motifEnrDT[[motifcolumname]]=="")
  logos <- paste("http://59.79.233.205/","/logos", "/","motifrow/", motifEnrDT[[motifcolumname]],".png", sep="")

  if(addHTML)
  {
    logos <- paste('<img src="', logos,
                   '" height="52" alt="',
                   motifEnrDT[[motifcolumname]], '"></img>', sep="")
  }
  logos[isNA] <- ""

  data.table::data.table(motifrowlogo=logos, motifEnrDT)
}


 #library(GenomicFeatures)
#library(GenomicAlignments)
#library(DESeq2)
#gtfFile <- file.path("/wrk/data/genome/yz_genome_data/mabamgenome/ma_genome.gff3")
#txdb <- makeTxDbFromGFF(file=gtfFile)
#isActiveSeq(txdb)
#gene_gr <- genes(txdb, columns="gene_id")
#bamfilename <- list.files("/wrk/zhu/work/FJ1_1__MN_bamboo/pre_processed/bamfiles/",pattern = ".bam$")
#bamfile <- paste0("/wrk/zhu/work/FJ1_1__MN_bamboo/pre_processed/bamfiles/",bamfilename)
#stopifnot(all(file.exists(bamfile)))
#samples_table <- read.delim("~/sampleinfo.csv")


##olap <- summarizeOverlaps(features = gene_gr, reads = bamLst, mode = "IntersectionStrict")

#colData(olap) <- DataFrame(samples_table)
#colnames(olap) <- samples_table$samplename


#deseq <- DESeqDataSet(olap, design = ~conduct)
#dese <- DESeq(deseq,test ="LRT",fitType='local',reduced=~1)
##res <- results(dese)
#summary(res)
#head(assay(olap), 3)

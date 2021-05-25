### ============================================================================
### EpiXprSData Package
### ----------------------------------------------------------------------------
###

## Download files
library(TCGAbiolinks)
library(EpiXprS)
rootDir <- getwd()
list <- c('TCGA-PRAD', 'TCGA-BRCA', 'TCGA-COAD', 'TCGA-KIRP', 'TCGA-UCEC', 'TCGA-LUAD',
  'TCGA-HNSC', 'TCGA-KIRC', 'TCGA-ESCA', 'TCGA-BLCA', 'TCGA-CESC', 'TCGA-CHOL',
  'TCGA-DLBC', 'TCGA-GBM', 'TCGA-LAML', 'TCGA-KICH', 'TCGA-LGG', 'TCGA-LIHC',
  'TCGA-LUSC', 'TCGA-MESO', 'TCGA-OV', 'TCGA-PAAD', 'TCGA-PCPG', 'TCGA-READ',
  'TCGA-SARC', 'TCGA-SKCM', 'TCGA-STAD', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM',
  'TCGA-UCS', 'TCGA-UVM', 'TCGA-ACC')

for(i in seq(cancers)){
    dir.create(paste0(rootDir,'/',list[i]))
### Download Clinical Data
query <- GDCquery(project = list[i], data.category = "Clinical",file.type='xml')
GDCdownload(query)

data <- GDCprepare_clinic(query, clinical.info = "patient")
clinical <- cbind(data@colData$sample,data@colData$race,data@colData$age_at_diagnosis)
clinical <- data.frame(clinical)
colnames(clinical) <- c('bcr_patient_barcode','race_list','age_at_initial_pathologic_diagnosis')
clinical$race_list <- toupper(clinical$race_list)
clinical$age_at_initial_pathologic_diagnosis <- round(as.numeric(as.character(clinical$age_at_initial_pathologic_diagnosis))/365.25)

### Download RNA counts and annotations
query.exp.hg38 <- GDCquery(project = paste0(list[i]),
                           data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification",
                           workflow.type = "HTSeq - Counts",
                           experimental.strategy = "RNA-Seq")
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(query = query.exp.hg38,
                     save = TRUE,
                     save.filename = "exp.Rda")

RNA_anno <- data.frame(expdat@rowRanges)
RNA_anno$seqnames <- substr(RNA_anno$seqnames,4,6)
RNA_clin <- cbind(RNA_anno$ensembl_gene_id,RNA_anno$seqnames,RNA_anno$start,RNA_anno$end,RNA_anno$external_gene_name)
colnames(RNA_clin) <- c('ensembl_id','chromosome_name','start_position','end_position','hgnc_symbol')

RNA <- data.frame(data@assays$data@listData)
colnames(RNA) <- clinical$bcr_patient_barcode
rownames(RNA) <- RNA_clin$ensembl_gene_id

### Download Methylation data
query_meth.hg38 <- GDCquery(project= paste0(list[i]),
                            data.category = "Gene expression",
                            platform = "Illumina Human Methylation 450")
GDCdownload(query_meth.hg38)
data.hg38 <- GDCprepare(query_meth.hg38,
                        save = TRUE,
                        save.filename = 'methy.Rda')

methy_anno <- data.frame(data.hg38@rowRanges)
methy_clin <- cbind(methy_anno$seqnames,methy_anno$start,methy_anno$Composite.Element.REF,methy_anno$Gene_Symbol)
methy_clin <- data.frame(methy_clin)
colnames(methy_clin) <-c('seqnames','start','Composite.Element.REF','Gene_Symbol')
methy <- data.frame(data.hg38@assays$data@listData)

r <- Construct(x = methy, y = RNA,
          clinical = clinical, method = 'Interaction', beta = FALSE,
          impute = FALSE)

tmp <- do.call(rbind,r[,2])
names(r) <- tmp$ensembl_gene
out <- r[,1]

for(i in seq(length(out))){
    out[[i]][[2]] <- tmp[i,]
    rownames(out[[i]][[1]]) <- out[[i]][[1]][,1]
    out[[i]][[1]] <- out[[i]][[1]][,-c(1)]
    out[[i]][[1]] <- data.frame(Weights = as.numeric(out[[i]][[1]]), row.names = names(out[[i]][[1]]))

}


paste0(list[i]) %>% save(file = paste0(list[i],'.rda'))

setwd(rootDir)
}

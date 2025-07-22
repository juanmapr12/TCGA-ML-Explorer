if(!exists("gdcGetURL_new", mode="function")){
  source("scripts/preparacion_entorno.R")
}


# Paquetes que nos van a hacer falta:
library(limma)
library(impute)
library(EnhancedVolcano)


# Estudiaremos aquí los datos de miRNAs
isoform_miRNAs_dir <- paste('metadatos', project, 'miRNAs', sep='/')

# Download miRNA data 
gdcRNADownload(project.id     = project, 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = isoform_miRNAs_dir)

# Query metadata from GDC graph
# Metadata associated with the isoform expression quantification file
metaMatrix.MIR <- gdcParseMetadata(project.id = project,
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

# Filtrar duplicados
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

# Filtrar las muestras non-Primary Tumor y non-Solid Tissue Normal en 
# los metadatos
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)

dim(metaMatrix.MIR)
colnames(metaMatrix.MIR)
str(metaMatrix.MIR)
head(metaMatrix.MIR)

# Define the function to merge all isoform quantification files into one df
merge_MIR <- function(metadata, fdir){
  filelist <- list.files(fdir, pattern="*.txt$", 
                         recursive = TRUE, full.names=TRUE)
  for (i in 1:length(filelist)){
    iname <- basename(filelist[i])
    isamplename <- metadata[metadata$file_name==iname, "sample"]
    if(length(isamplename) == 0)
      next
    idf <- read.csv(filelist[i], sep="\t", header=TRUE)
    idf <- idf %>% filter(startsWith(idf$miRNA_region,"mature"), idf$cross.mapped == "N")
    subset_idf <- idf[, c("miRNA_ID", "reads_per_million_miRNA_mapped")]
    subset_idf <- subset_idf %>% group_by(miRNA_ID) %>% summarise(conteo_normalizado = sum(reads_per_million_miRNA_mapped, na.rm = TRUE))
    names(subset_idf)[2] <- isamplename
    if (i==1){
      combined_df <- subset_idf
      rm(subset_idf)
    } else {
      combined_df <- merge(combined_df, subset_idf, by='miRNA_ID', all=TRUE)
      rm(subset_idf)
    }
  }
  rownames(combined_df) <- combined_df$miRNA_ID
  combined_df <- combined_df[,-which(names(combined_df) %in% "miRNA_ID")]
  return(combined_df)
}

MIRCounts <- merge_MIR(metaMatrix.MIR, isoform_miRNAs_dir)

# Pequeño análisis del conjunto de datos:
MIRCounts[1:5,1:20]
head(MIRCounts)
dim(MIRCounts)

MIRCounts_01 <- MIRCounts[grepl("01$", colnames(MIRCounts))]
MIRCounts_11 <- MIRCounts[grepl("11$", colnames(MIRCounts))]
dim(MIRCounts_01)
dim(MIRCounts_11)
unique(colnames(MIRCounts))

# (Para HNSC): 523 muestras cancerígenas y 44 muestras sanas


#### Puesta en marcha del conjunto de filtros ####

filtro_cerca_del_cero <- (rowSums(MIRCounts < 0.5, na.rm = TRUE) 
                         / ncol(MIRCounts)) < 0.5

valores_na_por_fila <- apply(MIRCounts, 1, function(x) sum(is.na(x)))
filtro_na <- valores_na_por_fila / ncol(MIRCounts) < 0.4

varianza <- apply(MIRCounts, 1, function(x) var(x, na.rm=TRUE))
filtro_desviacion <- varianza >= 0.5


filtros <- filtro_cerca_del_cero & filtro_na & filtro_desviacion
mirnas_tras_filtros <- MIRCounts[filtros,]


# Vamos a imputar los valores nulos mediante knn.

expr_mirnas_traspuesto <- t(mirnas_tras_filtros)
mirnas_tras_filtros <- t(impute.knn(as.matrix(expr_mirnas_traspuesto), 
                             k = 5)$data)

#### Puesta en marcha del análisis aplicado a los miRNAs ####

source("scripts/analisis_expr_dif.R")

# Busco ahora filtrar los genes más significativos (tanto en el
# logFC como en el p-valor) de ambos resultados:
resultados_global_significativos <- resultados_global[
  abs(resultados_global$logFC) >= 0.5 
  & resultados_global$adj.P.Val < 0.05,
]
resultados_pareadas_significativos <- resultados_pareadas[
  abs(resultados_pareadas$logFC) >= 0.5 
  & resultados_pareadas$adj.P.Val < 0.05,
]


# Y ahora simplemente consideramos la unión:
union_mirnas <- union(rownames(resultados_global_significativos),
                      rownames(resultados_pareadas_significativos))

mirnas_tras_filtros_y_dea <- 
  mirnas_tras_filtros[union_mirnas,]

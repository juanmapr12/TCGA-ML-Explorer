setwd("C:/Users/Usuario/OneDrive/Documentos/Archivos_R")
options(warn=-1)

project <- 'TCGA-HNSC'

# Descomentar estas líneas si no se tiene instalado el paquete 'BiocManager'
# ni la herramienta 'GDCRNATools:
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install('GDCRNATools')
  BiocManager::install('impute')
  BiocManager::install('EnhancedVolcano')
}
library(GDCRNATools)
library(moments)

# Librerías de tidyverse que nos van a servir:
library(ggplot2)
library(dplyr)
library(tidyr)

system("wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip")
unzip("gdc-client_v1.6.1_Ubuntu_x64.zip", exdir=".")

gdcGetURL_new <- function(project.id, data.type) {
  urlAPI <- 'https://api.gdc.cancer.gov/files?'
  if (data.type=='RNAseq') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'Gene Expression Quantification'
    workflow.type <- 'STAR - Counts'
  } else if (data.type=='miRNAs') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'Isoform Expression Quantification'
    workflow.type <- 'BCGSC miRNA Profiling'
  } else if (data.type=='Clinical') {
    data.category <- 'Clinical'
    data.type <- 'Clinical Supplement'
    workflow.type <- NA
  } else if (data.type=='pre-miRNAs') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'miRNA Expression Quantification'
    workflow.type <- 'BCGSC miRNA Profiling'
  }
  
  project <- paste('{"op":"in","content":{"field":"cases.',
                   'project.project_id","value":["', 
                   project.id, '"]}}', sep='')
  dataCategory <- paste('{"op":"in","content":{"field":"files.', 
                        'data_category","value":"', data.category, '"}}', sep='')
  dataType <- paste('{"op":"in","content":{"field":"files.data_type",',
                    '"value":"', data.type, '"}}', sep='')
  workflowType <- paste('{"op":"in","content":{"field":"files.',
                        'analysis.workflow_type","value":"', workflow.type, '"}}', sep='')
  
  
  if (is.na(workflow.type)) {
    dataFormat <- paste('{"op":"in","content":{"field":"files.',
                        'data_format","value":"', 'BCR XML', '"}}', sep='')
    content <- paste(project, dataCategory, dataType, dataFormat, sep=',')
  } else {
    content <- paste(project, dataCategory, dataType, 
                     workflowType, sep=',')
  }
  
  filters <- paste('filters=',URLencode(paste('{"op":"and","content":[', 
                                              content, ']}', sep='')),sep='')
  
  expand <- paste('analysis', 'analysis.input_files', 'associated_entities',
                  'cases', 'cases.diagnoses','cases.diagnoses.treatments', 
                  'cases.demographic', 'cases.project', 'cases.samples', 
                  'cases.samples.portions', 'cases.samples.portions.analytes', 
                  'cases.samples.portions.analytes.aliquots',
                  'cases.samples.portions.slides', sep=',')
  
  expand <- paste('expand=', expand, sep='')
  
  payload <- paste(filters, 'pretty=true', 'format=JSON', 
                   'size=10000', expand, sep='&')
  url <- paste(urlAPI, payload, sep='')
  
  return (url)
}

toolenv <- environment(get("gdcGetURL", envir = asNamespace("GDCRNATools")))
# Sacamos el entorno donde está definida la función 'gdcGetURL' dentro de el
# espacio de nombres de paquete 'GDCmiRNATools'.
# Dentro del entorno de la función puede haber variables y otras funciones que
# hagan referencia a ésta, de ahí que lo necesitemos.
# ¿Qué es un espacio de nombres de paquete? Es un ente que organiza y gestiona
# las funciones, variables y objetos dentro de un paquete.

unlockBinding("gdcGetURL", toolenv)
# Desbloqueamos variable que está bloqueada (por defecto) en el entorno con 
# el fin de modificarla

assignInNamespace("gdcGetURL", gdcGetURL_new, ns="GDCRNATools", envir=toolenv)
# Asignamos en el espacio de nombres o ns 'GDCRNATools' la nueva función
# definida sustituyendo a la antigua. Todo esto dentro del entorno.

assign("gdcGetURL", gdcGetURL_new)
# Asignación global (fuera de todo namespace)

lockBinding("gdcGetURL", toolenv)
# La volvemos a bloquear una vez realizado el cambio.

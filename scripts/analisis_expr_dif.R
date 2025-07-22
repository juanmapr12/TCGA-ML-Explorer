library(EnhancedVolcano)


# Análisis de expresión diferencial:

# Primero para el conjunto total de muestras:
muestras <- colnames(mirnas_tras_filtros)

analisis_expr_dif <- function(muestras, df){
  
  grupo <- vector(mode="character")
  
  for (i in 1:length(muestras)){
    tipo_muestra <- sapply(strsplit(muestras[i], "-"),
                           function(x) x[4])
    if (tipo_muestra == "01"){
      grupo <- c(grupo, "enfermo")
    }else{
      grupo <- c(grupo, "sano")
    }
  }
  
  grupo <- factor(grupo)
  grupo <- relevel(grupo, ref = "enfermo")
  design <- model.matrix(~ grupo)
  
  log_expr <- log2(df[,muestras] + 1)
  
  fit <- lmFit(log_expr, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf)
  
  return(res)
}


volcano <- function(res){
  print(EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'logFC',
                  y = 'adj.P.Val',
                  pCutoff = 0.05,
                  FCcutoff = 0.5))
}


resultados_global <- analisis_expr_dif(muestras, mirnas_tras_filtros)
volcano(resultados_global)



# Ahora lo mismo pero para muestras pareadas:

# Primer objetivo: sacar las muestras pareadas. Para ello, seguimos los
# siguientes pasos.

saca_pareadas <- function(){
  
  # 1.
  id_paciente <- sapply(strsplit(muestras, "-"), 
                        function(x) paste(x[1:3], collapse = "-"))
  
  # 2. 
  tipo_muestra <- sapply(strsplit(muestras, "-"), 
                         function(x) x[4])
  
  # 3.
  df <- data.frame(
    paciente = id_paciente,
    tipo = tipo_muestra
  )
  
  # 4. 
  muestras_pareadas <- df %>%
    group_by(paciente) %>%
    filter(all(c("01", "11") %in% tipo)) %>%
    ungroup() %>%
    arrange(paciente) %>%
    mutate(union = paste(paciente,tipo,sep="-")) %>%
    pull(union)
  
  return(muestras_pareadas)
}

muestras_pareadas <- saca_pareadas()

resultados_pareadas <- analisis_expr_dif(muestras_pareadas, mirnas_tras_filtros)
volcano(resultados_pareadas)



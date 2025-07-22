source("scripts/preprocesado_clinico.R")

# Librerías que vamos a necesitar:
library(randomForest)
library(randomForestSRC)
library(survival)
library(caret)

# Separamos en entrenamiento y prueba:
trainIndex <- createDataPartition(new_clinical_df$vital_status, p = 0.8, list = FALSE)
train_data <- new_clinical_df[trainIndex, ]
test_data <- new_clinical_df[-trainIndex, ]


# Modelo Random Forest

attach(train_data)
modelo_rf <- randomForest(vital_status ~ ., data = train_data, importance = TRUE)
print(modelo_rf)

# Evaluación en test
rf_preds <- predict(modelo_rf, newdata = test_data)
confusionMatrix(rf_preds, test_data$vital_status)

# Importancia de variables
varImpPlot(modelo_rf)



# Modelo Random Survival Forest

# #train_data_2 <- train_data %>% select(-time)
# rsf_model <- rfsrc(Surv(vital_status,time) ~ ., data = train_data, ntree = 1000, nodesize = 5, nsplit = 50, importance = TRUE)
# plot(rsf_model)
# print(rsf_model)

# RALLADA MENTAL. Aunque creo que nos puede venir bien para
# únicamente los miRNAs, no sé. Pero no todos obvio, quizás los que
# sobrevivan del análisis de supervivencia.

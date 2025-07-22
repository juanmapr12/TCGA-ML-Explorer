source("scripts/preprocesado_clinico.R")

# Librerías que vamos a necesitar:
library(randomForest)
library(randomForestSRC)
library(survival)
library(caret)

# Separamos en entrenamiento y prueba:
trainIndex <- createDataPartition(clinical_df$vital_status, p = 0.8, list = FALSE)
train_data <- clinical_df[trainIndex, ]
test_data <- clinical_df[-trainIndex, ]


# Modelo Random Forest

attach(train_data)
modelo_rf <- randomForest(vital_status ~ ., data = train_data, importance = TRUE)
print(modelo_rf)

# Evaluación en test
rf_preds <- predict(modelo_rf, newdata = test_data)
confusionMatrix(rf_preds, test_data$vital_status)

# Importancia de variables
varImpPlot(modelo_rf)

# Arbol elemental con libreria  rpart
# Curva ROC de la prediccion de un Arbol,  y AUC  Area Under Curve

# cargo las librerias que necesito
require("data.table")
require("rpart")
require("ROCR")

# Aqui se debe poner la carpeta de la materia de SU computadora local
setwd("~/buckets/b1/") # Establezco el Working Directory

# cargo el dataset
dataset <- fread("./datasets/competencia_2024.csv.gz")


# La curva ROC necesita clase binaria !
dataset[, clase_binaria1 := ifelse(clase_ternaria=="BAJA+2", "POS", "NEG" ) ]
dataset[, clase_ternaria := NULL ]

dtrain <- dataset[foto_mes == 202105] # defino donde voy a entrenar
dapply <- dataset[foto_mes == 202107] # defino donde voy a aplicar el modelo que tiene clase


# genero el modelo,  aqui se construye el arbol
modelo <- rpart(
        formula = "clase_binaria1 ~ .",
        data = dtrain, # los datos donde voy a entrenar
        xval = 0,
        cp =  -1, # esto significa no limitar la complejidad de los splits
        minsplit = 1078, # minima cantidad de registros para que se haga el split
        minbucket = 536, # tamaño minimo de una hoja
        maxdepth = 6
) # profundidad maxima del arbol



# aplico el modelo a los datos nuevos
prediccion <- predict(
        object = modelo,
        newdata = dapply,
        type = "prob"
)

curva_roc <-  ROCR::prediction(prediccion[, "POS"],
                 dapply$clase_binaria1,
                 label.ordering = c("NEG", "POS"))


dir.create("./exp/")
dir.create("./exp/ROC7001")


pdf( "./exp/ROC7001/ROC_rpart.pdf" )
perf <- performance(curva_roc,"tpr","fpr")
plot(perf)
dev.off()


perf <-  ROCR::performance(curva_roc, "auc")
AUC <- unlist(perf@y.values)

cat( "AUC = ",  AUC, "\n" )

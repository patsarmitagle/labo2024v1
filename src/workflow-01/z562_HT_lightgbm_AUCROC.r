# Experimentos Colaborativos Default
# Hyperparameter Tuning  lightgbm

# pensado para datasets con UNDERSAPLING de la clase mayoritaria

# limpio la memoria
rm(list = ls(all.names = TRUE)) # remove all objects
gc(full = TRUE) # garbage collection

require("data.table")
require("rlist")
require("yaml")

require("lightgbm")
require("ROCR")

# paquetes necesarios para la Bayesian Optimization
require("DiceKriging")
require("mlrMBO")

#------------------------------------------------------------------------------

options(error = function() {
  traceback(20)
  options(error = NULL)
  
  t <- format(Sys.time(), "%Y%m%d %H%M%S")
  cat( t, "\n",
    file = "z-Rabort.txt",
    append = TRUE
  )

  cat( t, "\n",
    file = "z-Rabort-hist.txt",
    append = TRUE
  )

  stop("exiting after script error")
})
#------------------------------------------------------------------------------

# Parametros del script
PARAM <- read_yaml( "parametros.yml" )

OUTPUT <- list()
#------------------------------------------------------------------------------

GrabarOutput <- function() {
  write_yaml(OUTPUT, file = "output.yml") # grabo OUTPUT
}
#------------------------------------------------------------------------------
# graba a un archivo los componentes de lista
# para el primer registro, escribe antes los titulos

exp_log <- function(
    reg, arch = NA, folder = "./exp/",
    ext = ".txt", verbose = TRUE) {
  archivo <- arch
  if (is.na(arch)) archivo <- paste0(folder, substitute(reg), ext)

  if (!file.exists(archivo)) # Escribo los titulos
    {
      linea <- paste0(
        "fecha\t",
        paste(list.names(reg), collapse = "\t"), "\n"
      )

      cat(linea, file = archivo)
    }

  linea <- paste0(
    format(Sys.time(), "%Y%m%d %H%M%S"), "\t", # la fecha y hora
    gsub(", ", "\t", toString(reg)), "\n"
  )

  cat(linea, file = archivo, append = TRUE) # grabo al archivo

  if (verbose) cat(linea) # imprimo por pantalla
}

#------------------------------------------------------------------------------
GLOBAL_arbol <- 0L
GLOBAL_gan_max <- -Inf
vcant_optima <- c()

#------------------------------------------------------------------------------

EstimarGanancia_lightgbm <- function(x) {
  gc()
  GLOBAL_iteracion <<- GLOBAL_iteracion + 1L
  OUTPUT$BO$iteracion_actual <<- GLOBAL_iteracion
  GrabarOutput()

  # hago la union de los parametros basicos y los moviles que vienen en x
  param_completo <- c(PARAM$lgb_basicos, x)

  param_completo$early_stopping_rounds <-
    as.integer(400 + 4 / param_completo$learning_rate)

  GLOBAL_arbol <<- 0L
  GLOBAL_gan_max <<- -Inf
  vcant_optima <<- c()
  set.seed(PARAM$lgb_semilla, kind = "L'Ecuyer-CMRG")
  modelo_train <- lgb.train(
    data = dtrain,
    valids = list(valid = dvalidate),
    eval = "auc",
    param = param_completo,
    verbose = -100
  )

  cat("\n")

  cant_corte <- vcant_optima[modelo_train$best_iter]

  # aplico el modelo a testing y calculo la ganancia
  prediccion <- predict(
    modelo_train,
    data.matrix(dataset_test[, campos_buenos, with = FALSE])
  )

  curva_roc <-  ROCR::prediction(prediccion,
                   dataset_test$clase01,
                   label.ordering = c(0, 1))

  perf <-  ROCR::performance(curva_roc, "auc")
  AUC_test <- unlist(perf@y.values)


  # logueo final
  ds <- list("cols" = ncol(dtrain), "rows" = nrow(dtrain))
  xx <- c(ds, copy(param_completo))

  xx$early_stopping_rounds <- NULL
  xx$num_iterations <- modelo_train$best_iter
  xx$AUC <- AUC_test
  xx$ganancia <- AUC_test
  xx$iteracion_bayesiana <- GLOBAL_iteracion

  exp_log(xx, arch = "BO_log.txt")

  # voy grabando las mejores column importance
  if (AUC_test > GLOBAL_AUC) {
    GLOBAL_AUC <<- AUC_test
    tb_importancia <- as.data.table(lgb.importance(modelo_train))

    fwrite(tb_importancia,
      file = paste0("impo_", sprintf("%03d", GLOBAL_iteracion), ".txt"),
      sep = "\t"
    )

    rm(tb_importancia)
    OUTPUT$BO$mejor$iteracion <<- GLOBAL_iteracion
    OUTPUT$BO$mejor$ganancia <<- GLOBAL_AUC
    OUTPUT$BO$mejor$arboles <<- modelo_train$best_iter
    GrabarOutput()
    exp_log(xx, arch = "BO_log_mejor.txt")
  }


  set.seed(PARAM$lgb_semilla, kind = "L'Ecuyer-CMRG")
  return(AUC_test)
}
#------------------------------------------------------------------------------

vcant_optima <- c()

EstimarGanancia_lightgbmCV <- function(x) {
  gc()
  GLOBAL_iteracion <<- GLOBAL_iteracion + 1L
  OUTPUT$BO$iteracion_actual <<- GLOBAL_iteracion
  GrabarOutput()

  param_completo <- c(PARAM$lgb_basicos, x)

  param_completo$early_stopping_rounds <-
    as.integer(400 + 4 / param_completo$learning_rate)

  vcant_optima <<- c()
  GLOBAL_arbol <<- 0L
  GLOBAL_gan_max <<- -Inf
  set.seed(PARAM$lgb_semilla, kind = "L'Ecuyer-CMRG")
  modelocv <- lgb.cv(
    data = dtrain,
    metric = "auc",
    param = param_completo,
    stratified = TRUE, # sobre el cross validation
    nfold = PARAM$lgb_crossvalidation_folds,
    verbose = -100
  )

  cat("\n")

  AUC_validation <- unlist(modelocv$record_evals$valid$auc$eval)[modelocv$best_iter]


  if (ktest == TRUE) {
    # debo recrear el modelo
    param_completo$early_stopping_rounds <- NULL
    param_completo$num_iterations <- modelocv$best_iter

    modelo <- lgb.train(
      data = dtrain,
      param = param_completo,
      verbose = -100
    )

    # aplico el modelo a testing y calculo la ganancia
    prediccion <- predict(
      modelo,
      data.matrix(dataset_test[, campos_buenos, with = FALSE])
    )

    curva_roc <-  ROCR::prediction(prediccion,
                     dataset_test$clase01,
                     label.ordering = c(0, 1))

    perf <-  ROCR::performance(curva_roc, "auc")
    AUC_test <- unlist(perf@y.values)

    rm(tbl)
    gc()
  }



  # voy grabando las mejores column importance
  if (AUC_test > GLOBAL_AUC) {
    GLOBAL_AUC <<- AUC_test

    param_impo <- copy(param_completo)
    param_impo$early_stopping_rounds <- 0
    param_impo$num_iterations <- modelocv$best_iter

    modelo <- lgb.train(
      data = dtrain,
      param = param_impo,
      verbose = -100
    )

    tb_importancia <- as.data.table(lgb.importance(modelo))

    fwrite(tb_importancia,
      file = paste0("impo_", GLOBAL_iteracion, ".txt"),
      sep = "\t"
    )

    rm(tb_importancia)

    OUTPUT$BO$mejor$AUC <<- GLOBAL_AUC
    OUTPUT$BO$mejor$ganancia <<- GLOBAL_AUC
    OUTPUT$BO$mejor$arboles <<- modelocv$best_iter
    GrabarOutput()
  }


  # logueo final
  ds <- list("cols" = ncol(dtrain), "rows" = nrow(dtrain))
  xx <- c(ds, copy(param_completo))

  xx$early_stopping_rounds <- NULL
  xx$num_iterations <- modelocv$best_iter
  xx$AUC <- AUC_test
  xx$ganancia <- AUC_test
  xx$iteracion_bayesiana <- GLOBAL_iteracion

  exp_log(xx, arch = "BO_log.txt")
  set.seed(PARAM$lgb_semilla, kind = "L'Ecuyer-CMRG")

  return(AUC_test)
}
#------------------------------------------------------------------------------

parametrizar  <- function( lparam )
{
  param_fijos  <- copy( lparam )
  hs  <- list()

  for( param  in  names( lparam ) )
  {
    if( length( lparam[[ param ]] ) > 1 )
    {
      desde  <- as.numeric( lparam[[ param ]][[1]]  )
      hasta  <- as.numeric( lparam[[ param ]][[2]]  )

      if( length( lparam[[ param ]] ) == 2 )
      {
         hs  <- append( hs,  
                        list( makeNumericParam( param, lower= desde, upper= hasta)  ) )
      } else {
         hs  <- append( hs, 
                        list( makeIntegerParam( param, lower= desde, upper= hasta) ) )
      }

      param_fijos[[ param ]] <- NULL  #lo quito 
    }
  }

  return( list( "param_fijos" =  param_fijos,
                "paramSet"    =  hs ) )
}

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Aqui empieza el programa
OUTPUT$PARAM <- PARAM
OUTPUT$time$start <- format(Sys.time(), "%Y%m%d %H%M%S")

# cargo las semillas
PARAM$lgb_semilla <- PARAM$semilla


# apertura de los parametros de LightGBM
#  en los que van fijos directo al LightGBM
#  y los que pasan a formar parte de la Bayesian Optimization
apertura  <- parametrizar( PARAM$lgb_param )
PARAM$lgb_basicos <- apertura$param_fijos
PARAM$bo_lgb <- makeParamSet( params= apertura$paramSet )

PARAM$lgb_basicos$seed <- PARAM$semilla

# cargo el dataset donde voy a entrenar
# esta en la carpeta del exp_input y siempre se llama  dataset.csv.gz
# cargo el dataset
PARAM$dataset <- paste0( "./", PARAM$input, "/dataset_training.csv.gz" )

dataset <- fread(PARAM$dataset)

dataset[, azar := NULL]

# Verificaciones
if (!("fold_train" %in% colnames(dataset))) {
  stop("Error, el dataset no tiene el campo fold_train \n")
}

if (!("fold_validate" %in% colnames(dataset))) {
  stop("Error, el dataset no tiene el campo fold_validate \n")
}

if (!("fold_test" %in% colnames(dataset))) {
  stop("Error, el dataset no tiene el campo fold_test  \n")
}

if (dataset[fold_train == 1, .N] == 0) {
  stop("Error, en el dataset no hay registros con fold_train==1 \n")
}


GrabarOutput()

cat(PARAM$exp_input,
  file = "TrainingStrategy.txt",
  append = FALSE
)

# defino la clase binaria clase01
dataset[, clase01 := ifelse(clase_ternaria == "CONTINUA", 0L, 1L)]


# los campos que se pueden utilizar para la prediccion
campos_buenos <- setdiff(
  copy(colnames(dataset)),
  c("clase01", "clase_ternaria", "fold_train", "fold_validate", "fold_test")
)

# la particion de train siempre va
dtrain <- lgb.Dataset(
  data = data.matrix(dataset[fold_train == 1, campos_buenos, with = FALSE]),
  label = dataset[fold_train == 1, clase01],
  weight = dataset[
    fold_train == 1,
    ifelse(clase_ternaria == "BAJA+2", 1.0000001,
      ifelse(clase_ternaria == "BAJA+1", 1.0, 1.0)
    )
  ],
  free_raw_data = FALSE
)

OUTPUT$train$ncol <- ncol(dtrain)
OUTPUT$train$nrow <- nrow(dtrain)
OUTPUT$train$periodos <- dataset[fold_train == 1, length(unique(foto_mes))]

kvalidate <- FALSE
ktest <- FALSE
kcrossvalidation <- TRUE

# Si hay que hacer validacion
if (dataset[fold_train == 0 & fold_test == 0 & fold_validate == 1, .N] > 0) {
  kcrossvalidation <- FALSE
  kvalidate <- TRUE
  dvalidate <- lgb.Dataset(
    data = data.matrix(dataset[fold_validate == 1, campos_buenos, with = FALSE]),
    label = dataset[fold_validate == 1, clase01],
    weight = dataset[
      fold_validate == 1,
      ifelse(clase_ternaria == "BAJA+2", 1.0000001,
        ifelse(clase_ternaria == "BAJA+1", 1.0, 1.0)
      )
    ],
    free_raw_data = FALSE
  )

  OUTPUT$validate$ncol <- ncol(dvalidate)
  OUTPUT$validate$nrow <- nrow(dvalidate)

  OUTPUT$validate$periodos <- dataset[
    fold_validate == 1,
    length(unique(foto_mes))
  ]
}


# Si hay que hacer testing
if (dataset[fold_train == 0 & fold_validate == 0 & fold_test == 1, .N] > 0) {
  ktest <- TRUE
  campos_buenos_test <- setdiff(
    copy(colnames(dataset)),
    c("fold_train", "fold_validate", "fold_test")
  )

  dataset_test <- dataset[fold_test == 1, campos_buenos_test, with = FALSE]

  OUTPUT$test$ncol <- ncol(dataset_test)
  OUTPUT$test$nrow <- nrow(dataset_test)
  OUTPUT$test$periodos <- dataset_test[, length(unique(foto_mes))]
}



rm(dataset)
gc()


# si ya existe el archivo log, traigo hasta donde procese
if (file.exists("BO_log.txt")) {
  tabla_log <- fread("BO_log.txt")
  GLOBAL_iteracion <- nrow(tabla_log)
  GLOBAL_AUC <- tabla_log[, max(ganancia)]
  rm(tabla_log)
} else {
  GLOBAL_iteracion <- 0L
  GLOBAL_AUC <- -Inf
}


# Aqui comienza la configuracion de mlrMBO


OUTPUT$crossvalidation <- kcrossvalidation
GrabarOutput()

# deobo hacer cross validation o  Train/Validate/Test
if (kcrossvalidation) {
  funcion_optimizar <- EstimarGanancia_lightgbmCV
} else {
  funcion_optimizar <- EstimarGanancia_lightgbm
}


configureMlr(show.learner.output = FALSE)

# configuro la busqueda bayesiana,  los hiperparametros que se van a optimizar
# por favor, no desesperarse por lo complejo
obj.fun <- makeSingleObjectiveFunction(
  fn = funcion_optimizar, # la funcion que voy a maximizar
  minimize = FALSE, # estoy Maximizando la ganancia
  noisy = TRUE,
  par.set = PARAM$bo_lgb, # definido al comienzo del programa
  has.simple.signature = FALSE # paso los parametros en una lista
)

# archivo donde se graba y cada cuantos segundos
ctrl <- makeMBOControl(
  save.on.disk.at.time = 600,
  save.file.path = "bayesiana.RDATA"
)

ctrl <- setMBOControlTermination(ctrl,
  iters = PARAM$bo_iteraciones
) # cantidad de iteraciones

ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())

# establezco la funcion que busca el maximo
surr.km <- makeLearner("regr.km",
  predict.type = "se",
  covtype = "matern3_2",
  control = list(trace = TRUE)
)

surr.km <- makeLearner("regr.km",
  predict.type = "se",
  covtype = "matern3_2",
  optim.method = "BFGS",
  nugget.estim = TRUE,
  jitter = TRUE,
  control = list(trace = TRUE)
)


# Aqui inicio la optimizacion bayesiana
set.seed(PARAM$lgb_semilla, kind = "L'Ecuyer-CMRG")
if (!file.exists("bayesiana.RDATA")) {
  run <- mbo(obj.fun, learner = surr.km, control = ctrl)
} else {
  # si ya existe el archivo RDATA,
  # debo continuar desde el punto hasta donde llegue
  #  usado para cuando se corta la virtual machine
  run <- mboContinue("bayesiana.RDATA") # retomo en caso que ya exista
}

#------------------------------------------------------------------------------
BO_log <- fread("BO_log.txt")
OUTPUT$ganancia_max <- BO_log[, max(ganancia, na.rm = TRUE)]

OUTPUT$time$end <- format(Sys.time(), "%Y%m%d %H%M%S")
GrabarOutput()

# dejo la marca final
cat(format(Sys.time(), "%Y%m%d %H%M%S"), "\n",
  file = "z-Rend.txt",
  append = TRUE
)

# Optimizacion Bayesiana de hiperparametros de  rpart
# 10-repeated 5-fold cross validation

# limpio la memoria
rm(list = ls()) # remove all objects
gc() # garbage collection

require("data.table")
require("rlist")
require("primes")

require("rpart")
require("parallel")

# paquetes necesarios para la Bayesian Optimization
require("DiceKriging")
require("mlrMBO")


# Defino la  Optimizacion Bayesiana

# cantidad de iteraciones de la Optimizacion Bayesiana
PARAM <- list()

PARAM$xval_repeated <- 10
PARAM$BO_iter <- 50 #cantidad de iteraciones de la Bayesian Optimization

# la letra L al final de 1L significa ENTERO
PARAM$hs <- makeParamSet(
    makeNumericParam("cp", lower = -1, upper = 0.1),
    makeIntegerParam("minsplit", lower = 1L, upper = 8000L),
    makeIntegerParam("minbucket", lower = 1L, upper = 4000L),
    makeIntegerParam("maxdepth", lower = 3L, upper = 20L),
    forbidden = quote(minbucket > 0.5 * minsplit)
)
# minbuket NO PUEDE ser mayor que la mitad de minsplit


#------------------------------------------------------------------------------
# genera un vector de numeros primos aleatorios
# partiendo de MI primer semilla

generarprimos <- function( qty ){

  # recupero MI primer semilla
  tabla_semillas <- fread( "./datasets//mis_semillas.txt" )
  mi_semilla <- tabla_semillas[ 1, semilla ]  # 1 es mi primer semilla

  primos_universo  <- generate_primes(min=100000, max=1000000)  #genero TODOS los numeros primos entre 100k y 1M
  set.seed( mi_semilla ) #seteo la semilla que controla al sample de los primos

  return( sample(primos_universo, qty) )
}

#------------------------------------------------------------------------------
# graba a un archivo los componentes de lista
# para el primer registro, escribe antes los titulos

loguear <- function(reg, arch = NA, folder = "./work/", ext = ".txt",
                    verbose = TRUE) {
  archivo <- arch
  if (is.na(arch)) archivo <- paste0(folder, substitute(reg), ext)

  # Escribo los titulos
  if (!file.exists(archivo)) {
    linea <- paste0(
      "fecha\t",
      paste(list.names(reg), collapse = "\t"), "\n"
    )

    cat(linea, file = archivo)
  }

  # la fecha y hora
  linea <- paste0(
    format(Sys.time(), "%Y%m%d %H%M%S"), "\t",
    gsub(", ", "\t", toString(reg)), "\n"
  )

  # grabo al archivo
  cat(linea, file = archivo, append = TRUE)

  # imprimo por pantalla
  if (verbose) cat(linea)
}
#------------------------------------------------------------------------------
# particionar agrega una columna llamada fold a un dataset
#   que consiste en una particion estratificada segun agrupa
# particionar( data=dataset, division=c(70,30),
#  agrupa=clase_ternaria, seed=semilla)   crea una particion 70, 30
# particionar( data=dataset, division=c(1,1,1,1,1),
#  agrupa=clase_ternaria, seed=semilla)   divide el dataset en 5 particiones

particionar <- function(data, division, agrupa = "", campo = "fold",
                        start = 1, seed = NA) {
  if (!is.na(seed)) set.seed(seed)

  bloque <- unlist(mapply(
    function(x, y) {
      rep(y, x)
    }, division,
    seq(from = start, length.out = length(division))
  ))

  data[, (campo) := sample(rep(bloque, ceiling(.N / length(bloque))))[1:.N],
    by = agrupa
  ]
}
#------------------------------------------------------------------------------
# fold_test  tiene el numero de fold que voy a usar para testear,
#  entreno en el resto de los folds
# param tiene los hiperparametros del arbol

ArbolSimple <- function(fold_test, data, param) {
  # genero el modelo
  # entreno en todo MENOS el fold_test que uso para testing
  modelo <- rpart("clase_ternaria ~ .",
    data = data[fold != fold_test, ],
    xval = 0,
    control = param
  )

  # aplico el modelo a los datos de testing
  # aplico el modelo sobre los datos de testing
  # quiero que me devuelva probabilidades
  prediccion <- predict(modelo,
    data[fold == fold_test, ],
    type = "prob"
  )

  # esta es la probabilidad de baja
  prob_baja2 <- prediccion[, "BAJA+2"]

  # calculo la ganancia
  ganancia_testing <- data[fold == fold_test][
    prob_baja2 > 1 / 40,
    sum(ifelse(clase_ternaria == "BAJA+2",
      117000, -3000
    ))
  ]

  # esta es la ganancia sobre el fold de testing, NO esta normalizada
  return(ganancia_testing)
}
#------------------------------------------------------------------------------

ArbolesCrossValidation <- function(data, param, qfolds, pagrupa, semilla) {
  # generalmente  c(1, 1, 1, 1, 1 )  cinco unos
  divi <- rep(1, qfolds)

  # particiono en dataset en folds
  particionar(data, divi, seed = semilla, agrupa = pagrupa)

  ganancias <- mcmapply(ArbolSimple,
    seq(qfolds), # 1 2 3 4 5
    MoreArgs = list(data, param),
    SIMPLIFY = FALSE,
    mc.cores = qfolds
  )

  data[, fold := NULL]

  # devuelvo la primer ganancia y el promedio
  # promedio las ganancias
  ganancia_promedio <- mean(unlist(ganancias))
  # aqui normalizo la ganancia
  ganancia_promedio_normalizada <- ganancia_promedio * qfolds

  return(ganancia_promedio_normalizada)
}
#------------------------------------------------------------------------------
# esta funcion solo puede recibir los parametros que se estan optimizando
# el resto de los parametros, lamentablemente se pasan como variables globales

EstimarGanancia <- function(x) {
  GLOBAL_iteracion <<- GLOBAL_iteracion + 1

  # creo un vector vacio donde ire agregando las ganancias
  vganancias <- c()

  # recorro las semillas
  for( semillita in kvector_semillas )
  {
    xval_folds <- 5
    # param= x los hiperparametros del arbol
    # qfolds= xval_folds  la cantidad de folds
    ganancia <- ArbolesCrossValidation(dataset,
      param = x,
      qfolds = xval_folds,
      pagrupa = "clase_ternaria",
      semilla = semillita
    )

    # agrego al vector
    vganancias <- c( vganancias, ganancia) 
  }

  ganancia_promedio <- mean( vganancias )
  # logueo
  xx <- x
  xx$xval_folds <- xval_folds
  xx$xval_repe <- length( kvector_semillas )
  xx$ganancia <- ganancia_promedio
  xx$iteracion <- GLOBAL_iteracion

  # si es ganancia_promedio superadora la almaceno en mejor
  if( ganancia_promedio > GLOBAL_mejor ) {
    GLOBAL_mejor <<- ganancia_promedio
    Sys.sleep(2)
    loguear(xx, arch = archivo_log_mejor)
  }
  

  Sys.sleep(2)
  loguear(xx, arch = archivo_log)

  return(ganancia_promedio)
}
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Aqui empieza el programa

# Establezco el Working Directory
setwd("~/buckets/b1/")

# genero PARAM$xval_repeated   numeros primos
kvector_semillas <- generarprimos( PARAM$xval_repeated )

# cargo los datos
dataset <- fread("./datasets/dataset_pequeno.csv")
# entreno en 202107
dataset <- dataset[foto_mes==202107]


# creo la carpeta donde va el experimento
#  HT  representa  Hiperparameter Tuning
dir.create("./exp/", showWarnings = FALSE)
dir.create("./exp/HT3330/", showWarnings = FALSE)

# Establezco el Working Directory DEL EXPERIMENTO
setwd("./exp/HT3330/")


archivo_log <- "HT333.txt"
archivo_log_mejor <- "HT333_mejor.txt"
archivo_BO <- "HT333.RDATA"

# leo si ya existe el log
#  para retomar en caso que se se corte el programa
GLOBAL_iteracion <- 0
GLOBAL_mejor <- -Inf

if (file.exists(archivo_log)) {
  tabla_log <- fread(archivo_log)
  GLOBAL_iteracion <- nrow(tabla_log)
  GLOBAL_mejor <- tabla_log[, max(ganancia)]
}



# Aqui comienza la configuracion de la Bayesian Optimization

funcion_optimizar <- EstimarGanancia

configureMlr(show.learner.output = FALSE)

# configuro la busqueda bayesiana,
#  los hiperparametros que se van a optimizar
# por favor, no desesperarse por lo complejo
# minimize= FALSE estoy Maximizando la ganancia
obj.fun <- makeSingleObjectiveFunction(
  fn = funcion_optimizar,
  minimize = FALSE,
  noisy = TRUE,
  par.set = PARAM$hs,
  has.simple.signature = FALSE
)

ctrl <- makeMBOControl(
  save.on.disk.at.time = 600,
  save.file.path = archivo_BO
)

ctrl <- setMBOControlTermination(ctrl, iters = PARAM$BO_iter)
ctrl <- setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())

surr.km <- makeLearner("regr.km",
  predict.type = "se",
  covtype = "matern3_2", control = list(trace = TRUE)
)

# inicio la optimizacion bayesiana
if (!file.exists(archivo_BO)) {
  run <- mbo(
    fun = obj.fun,
    learner = surr.km,
    control = ctrl
  )
} else {
  run <- mboContinue(archivo_BO)
}
# retomo en caso que ya exista

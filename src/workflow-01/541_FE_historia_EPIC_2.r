# Experimentos Colaborativos Default
# Workflow  Feature Engineering historico

# limpio la memoria
rm(list = ls(all.names = TRUE)) # remove all objects
gc(full = TRUE) # garbage collection

require("data.table")
require("yaml")
require("Rcpp")
library("stats")
require("ranger")
require("randomForest") # solo se usa para imputar nulos
library(zoo)

require("lightgbm")


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
# se calculan para los 6 meses previos el minimo, maximo y
#  tendencia calculada con cuadrados minimos
# la formula de calculo de la tendencia puede verse en
#  https://stats.libretexts.org/Bookshelves/Introductory_Statistics/Book%3A_Introductory_Statistics_(Shafer_and_Zhang)/10%3A_Correlation_and_Regression/10.04%3A_The_Least_Squares_Regression_Line
# para la maxíma velocidad esta funcion esta escrita en lenguaje C,
# y no en la porqueria de R o Python

cppFunction("NumericVector fhistC(NumericVector pcolumna, IntegerVector pdesde )
{
  /* Aqui se cargan los valores para la regresion */
  double  x[100] ;
  double  y[100] ;

  int n = pcolumna.size();
  NumericVector out( 7*n );

  for(int i = 0; i < n; i++)
  {
    //lag
    if( pdesde[i]-1 < i )  out[ i + 4*n ]  =  pcolumna[i-1] ;
    else                   out[ i + 4*n ]  =  NA_REAL ;


    int  libre    = 0 ;
    int  xvalor   = 1 ;

    for( int j= pdesde[i]-1;  j<=i; j++ )
    {
       double a = pcolumna[j] ;

       if( !R_IsNA( a ) )
       {
          y[ libre ]= a ;
          x[ libre ]= xvalor ;
          libre++ ;
       }

       xvalor++ ;
    }

    /* Si hay al menos dos valores */
    if( libre > 1 )
    {
      double  xsum  = x[0] ;
      double  ysum  = y[0] ;
      double  xysum = xsum * ysum ;
      double  xxsum = xsum * xsum ;
      double  vmin  = y[0] ;
      double  vmax  = y[0] ;

      for( int h=1; h<libre; h++)
      {
        xsum  += x[h] ;
        ysum  += y[h] ;
        xysum += x[h]*y[h] ;
        xxsum += x[h]*x[h] ;

        if( y[h] < vmin )  vmin = y[h] ;
        if( y[h] > vmax )  vmax = y[h] ;
      }

      out[ i ]  =  (libre*xysum - xsum*ysum)/(libre*xxsum -xsum*xsum) ;
      out[ i + n ]    =  vmin ;
      out[ i + 2*n ]  =  vmax ;
      out[ i + 3*n ]  =  ysum / libre ;      
      
      // Calculate autocorrelation for a chosen lag (k)
      // Autocorrelación
      double autocorr_sum = 0.0;
      for(int h = 0; h < libre - 1; h++)
      {
        autocorr_sum += (y[h] - out[i + 3 * n]) * (y[h + 1] - out[i + 3 * n]);
      }
      out[i + 4 * n] = autocorr_sum / (libre - 1);
      
      // Momento (varianza)
      double moment_sum = 0.0;
      for(int h = 0; h < libre; h++)
      {
        moment_sum += pow(y[h] - out[i + 3 * n], 2);
      }
      out[i + 5 * n] = moment_sum / (libre - 1);
    }
    else
    {
      out[ i       ]  =  NA_REAL ;
      out[ i + n   ]  =  NA_REAL ;
      out[ i + 2*n ]  =  NA_REAL ;
      out[ i + 3*n ]  =  NA_REAL ;
      out[i + 4 * n] = NA_REAL; // Autocorrelación
      out[i + 5 * n] = NA_REAL; // Momento
    }
  }

  return  out;
}")

## calcula todos los ratios entre campos monetarios
RatiosEpico <- function(
  dataset, cols
){
  for (i in 1:(length(cols) -1)){
    for (j in (i+1):length(cols)){
      campo_1 <- cols[i]
      campo_2 <- cols[j]

      nueva_col <- paste(campo_1, campo_2, "ratio", sep="_")
      ratio <- ifelse(dataset[[campo_2]] == 0,
                sign(dataset[[campo_1]]) * .Machine$double.xmax,
                dataset[[campo_1]] / dataset[[campo_2]])
      dataset[nueva_col] <- ratio
    }
  }
}
RatiosEpico <- function(dataset, cols) {
  
  for (i in 1:(length(cols) - 1)) {
    for (j in (i + 1):length(cols)) {
      campo_1 <- cols[i]
      campo_2 <- cols[j]

      nueva_col <- paste(campo_1, campo_2, "ratio", sep = "_")
      
      # Calcular el ratio y manejar divisiones por cero
      dataset[, (nueva_col) := ifelse(get(campo_2) == 0,
                                      sign(get(campo_1)) * .Machine$double.xmax,
                                      get(campo_1) / get(campo_2))]
    }
  }
  #return(dataset)
}

RatiosEpicoDiv0HaceNA <- function(
  dataset, cols
){
  for (i in 1:(length(cols) -1)){
    for (j in (i+1):length(cols)){
      campo_1 <- cols[i]
      campo_2 <- cols[j]

      nueva_col <- paste(campo_1, campo_2, "ratio", sep="_")
      ratio <- ifelse(dataset[[campo_2]] == 0,
                NaN,
                dataset[[campo_1]] / dataset[[campo_2]])
      dataset[nueva_col] <- ratio
    }
  }
}

#------------------------------------------------------------------------------
# calcula la tendencia de las variables cols de los ultimos 6 meses
# la tendencia es la pendiente de la recta que ajusta por cuadrados minimos
# La funcionalidad de ratioavg es autoria de  Daiana Sparta,  UAustral  2021

TendenciaYmuchomas <- function(
    dataset, cols, ventana = 6, tendencia = TRUE,
    minimo = TRUE, maximo = TRUE, promedio = TRUE,
    ratioavg = FALSE, ratiomax = FALSE, autocorr = FALSE) {
  gc()
  # Esta es la cantidad de meses que utilizo para la historia
  ventana_regresion <- ventana

  last <- nrow(dataset)

  # creo el vector_desde que indica cada ventana
  # de esta forma se acelera el procesamiento ya que lo hago una sola vez
  vector_ids <- dataset[ , get( PARAM$dataset_metadata$entity_id) ]

  vector_desde <- seq(
    -ventana_regresion + 2,
    nrow(dataset) - ventana_regresion + 1
  )

  vector_desde[1:ventana_regresion] <- 1

  for (i in 2:last) {
    if (vector_ids[i - 1] != vector_ids[i]) {
      vector_desde[i] <- i
    }
  }
  for (i in 2:last) {
    if (vector_desde[i] < vector_desde[i - 1]) {
      vector_desde[i] <- vector_desde[i - 1]
    }
  }

  for (campo in cols) {
    nueva_col <- fhistC(dataset[, get(campo)], vector_desde)

    if (autocorr) {
        dataset[, paste0(campo, "_autocorr", ventana) :=
          nueva_col[(4 * last + 1):(5 * last)]]
        dataset[, paste0(campo, "_momentum", ventana) :=
                nueva_col[(5 * last + 1):(6 * last)]]
    }

    if (tendencia) {
      dataset[, paste0(campo, "_tend", ventana) :=
        nueva_col[(0 * last + 1):(1 * last)]]
    }

    if (minimo) {
      dataset[, paste0(campo, "_min", ventana) :=
        nueva_col[(1 * last + 1):(2 * last)]]
    }

    if (maximo) {
      dataset[, paste0(campo, "_max", ventana) :=
        nueva_col[(2 * last + 1):(3 * last)]]
    }

    if (promedio) {
      dataset[, paste0(campo, "_avg", ventana) :=
        nueva_col[(3 * last + 1):(4 * last)]]
    }

    if (ratioavg) {
      dataset[, paste0(campo, "_ratioavg", ventana) :=
        get(campo) / nueva_col[(3 * last + 1):(4 * last)]]
    }

    if (ratiomax) {
      dataset[, paste0(campo, "_ratiomax", ventana) :=
        get(campo) / nueva_col[(2 * last + 1):(3 * last)]]
    }
  }
  return (dataset)
}
#------------------------------------------------------------------------------
# agrega al dataset nuevas variables {0,1}
#  que provienen de las hojas de un Random Forest

AgregaVarRandomForest <- function(
    num.trees, max.depth,
    min.node.size, mtry, semilla) {
  gc()
  dataset[, clase01 := ifelse(clase_ternaria == "CONTINUA", 0, 1)]

  campos_buenos <- setdiff(colnames(dataset), c("clase_ternaria"))

  dataset_rf <- copy(dataset[, campos_buenos, with = FALSE])
  set.seed(semilla, kind = "L'Ecuyer-CMRG")
  azar <- runif(nrow(dataset_rf))

  dataset_rf[, entrenamiento :=
    as.integer(foto_mes >= 202007 & foto_mes <= 202103 &
      (clase01 == 1 | azar < 0.10))]

  # imputo los nulos, ya que ranger no acepta nulos
  # Leo Breiman, ¿por que le temias a los nulos?
  set.seed(semilla, kind = "L'Ecuyer-CMRG")
  dataset_rf <- na.roughfix(dataset_rf)

  campos_buenos <- setdiff(
    colnames(dataset_rf),
    c("clase_ternaria", "entrenamiento")
  )

  set.seed(semilla, kind = "L'Ecuyer-CMRG")
  modelo <- ranger(
    formula = "clase01 ~ .",
    data = dataset_rf[entrenamiento == 1L, campos_buenos, with = FALSE],
    classification = TRUE,
    probability = FALSE,
    num.trees = num.trees,
    max.depth = max.depth,
    min.node.size = min.node.size,
    mtry = mtry,
    seed = semilla,
    num.threads = 1
  )

  rfhojas <- predict(
    object = modelo,
    data = dataset_rf[, campos_buenos, with = FALSE],
    predict.all = TRUE, # entrega la prediccion de cada arbol
    type = "terminalNodes" # entrega el numero de NODO el arbol
  )

  for (arbol in 1:num.trees) {
    hojas_arbol <- unique(rfhojas$predictions[, arbol])

    for (pos in 1:length(hojas_arbol)) {
      # el numero de nodo de la hoja, estan salteados
      nodo_id <- hojas_arbol[pos]
      dataset[, paste0(
        "rf_", sprintf("%03d", arbol),
        "_", sprintf("%03d", nodo_id)
      ) := 0L]

      dataset[
        which(rfhojas$predictions[, arbol] == nodo_id, ),
        paste0(
          "rf_", sprintf("%03d", arbol),
          "_", sprintf("%03d", nodo_id)
        ) := 1L
      ]
    }
  }

  rm(dataset_rf)
  dataset[, clase01 := NULL]

  gc()
}
#------------------------------------------------------------------------------
VPOS_CORTE <- c()

fganancia_lgbm_meseta <- function(probs, datos) {
  vlabels <- get_field(datos, "label")
  vpesos <- get_field(datos, "weight")

  tbl <- as.data.table(list(
    "prob" = probs,
    "gan" = ifelse(vlabels == 1 & vpesos > 1, 117000, -3000)
  ))

  setorder(tbl, -prob)
  tbl[, posicion := .I]
  tbl[, gan_acum := cumsum(gan)]
  setorder(tbl, -gan_acum) # voy por la meseta

  gan <- mean(tbl[1:500, gan_acum]) # meseta de tamaño 500

  pos_meseta <- tbl[1:500, median(posicion)]
  VPOS_CORTE <<- c(VPOS_CORTE, pos_meseta)

  return(list(
    "name" = "ganancia",
    "value" = gan,
    "higher_better" = TRUE
  ))
}
#------------------------------------------------------------------------------
# Elimina del dataset las variables que estan por debajo
#  de la capa geologica de canaritos
# se llama varias veces, luego de agregar muchas variables nuevas,
#  para ir reduciendo la cantidad de variables
# y así hacer lugar a nuevas variables importantes

GVEZ <- 1

CanaritosAsesinos <- function(
    canaritos_ratio = 0.2,
    canaritos_desvios = 3.0, canaritos_semilla = 999983) {
  gc()
  dataset[, clase01 := ifelse(clase_ternaria == "CONTINUA", 0, 1)]

  set.seed(canaritos_semilla, kind = "L'Ecuyer-CMRG")
  for (i in 1:(ncol(dataset) * canaritos_ratio)) {
    dataset[, paste0("canarito", i) := runif(nrow(dataset))]
  }

  campos_buenos <- setdiff(
    colnames(dataset),
    c( campitos, "clase01")
  )

  azar <- runif(nrow(dataset))

  dataset[, entrenamiento :=
    foto_mes >= 202009 & foto_mes <= 202103 & (clase01 == 1 | azar < 0.10)]

  dtrain <- lgb.Dataset(
    data = data.matrix(dataset[entrenamiento == TRUE, campos_buenos, with = FALSE]),
    label = dataset[entrenamiento == TRUE, clase01],
    weight = dataset[
      entrenamiento == TRUE,
      ifelse(clase_ternaria == "BAJA+2", 1.0000001, 1.0)
    ],
    free_raw_data = FALSE
  )

  dvalid <- lgb.Dataset(
    data = data.matrix(dataset[foto_mes == 202105, campos_buenos, with = FALSE]),
    label = dataset[foto_mes == 202105, clase01],
    weight = dataset[
      foto_mes == 202105,
      ifelse(clase_ternaria == "BAJA+2", 1.0000001, 1.0)
    ],
    free_raw_data = FALSE
  )


  param <- list(
    objective = "binary",
    metric = "custom",
    first_metric_only = TRUE,
    boost_from_average = TRUE,
    feature_pre_filter = FALSE,
    verbosity = -100,
    seed = canaritos_semilla,
    max_depth = -1, # -1 significa no limitar,  por ahora lo dejo fijo
    min_gain_to_split = 0.0, # por ahora, lo dejo fijo
    lambda_l1 = 0.0, # por ahora, lo dejo fijo
    lambda_l2 = 0.0, # por ahora, lo dejo fijo
    max_bin = 31, # por ahora, lo dejo fijo
    num_iterations = 9999, # un numero grande, lo limita early_stopping_rounds
    force_row_wise = TRUE, # para que los alumnos no se atemoricen con  warning
    learning_rate = 0.065,
    feature_fraction = 1.0, # lo seteo en 1
    min_data_in_leaf = 260,
    num_leaves = 60,
    early_stopping_rounds = 200,
    num_threads = 1
  )

  set.seed(canaritos_semilla, kind = "L'Ecuyer-CMRG")
  modelo <- lgb.train(
    data = dtrain,
    valids = list(valid = dvalid),
    eval = fganancia_lgbm_meseta,
    param = param,
    verbose = -100
  )

  tb_importancia <- lgb.importance(model = modelo)
  tb_importancia[, pos := .I]

  fwrite(tb_importancia,
    file = paste0("impo_", GVEZ, ".txt"),
    sep = "\t"
  )

  GVEZ <<- GVEZ + 1

  umbral <- tb_importancia[
    Feature %like% "canarito",
    median(pos) + canaritos_desvios * sd(pos)
  ] # Atencion corto en la mediana mas desvios!!

  col_utiles <- tb_importancia[
    pos < umbral & !(Feature %like% "canarito"),
    Feature
  ]

  col_utiles <- unique(c(
    col_utiles,
    c(campitos, "mes")
  ))

  col_inutiles <- setdiff(colnames(dataset), col_utiles)

  dataset[, (col_inutiles) := NULL]
}
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Aqui empieza el programa
OUTPUT$PARAM <- PARAM
OUTPUT$time$start <- format(Sys.time(), "%Y%m%d %H%M%S")

PARAM$RandomForest$semilla <- PARAM$semilla
PARAM$CanaritosAsesinos1$semilla <- PARAM$semilla
PARAM$CanaritosAsesinos2$semilla <- PARAM$semilla

# cargo el dataset donde voy a entrenar
# esta en la carpeta del exp_input y siempre se llama  dataset.csv.gz
# cargo el dataset
PARAM$dataset <- paste0( "./", PARAM$input, "/dataset.csv.gz" )
PARAM$dataset_metadata <- read_yaml( paste0( "./", PARAM$input, "/dataset_metadata.yml" ) )

dataset <- fread(PARAM$dataset)


colnames(dataset)[which(!(sapply(dataset, typeof) %in% c("integer", "double")))]


GrabarOutput()

#--------------------------------------
# estas son las columnas a las que se puede agregar
#  lags o media moviles ( todas menos las obvias )

campitos <- c( PARAM$dataset_metadata$primarykey,
  PARAM$dataset_metadata$entity_id,
  PARAM$dataset_metadata$periodo,
  PARAM$dataset_metadata$clase )

campitos <- unique( campitos )

cols_lagueables <- copy(setdiff(
  colnames(dataset),
  PARAM$dataset_metadata
))

# ordeno el dataset por <numero_de_cliente, foto_mes> para poder hacer lags
#  es MUY  importante esta linea
# ordeno dataset
setorderv(dataset, PARAM$dataset_metadata$primarykey)


cols_monetarios <- cols_lagueables
cols_monetarios <- cols_monetarios[cols_monetarios %like%
  "^(m|Visa_m|Master_m|vm_m|cliente_edad|cliente_antigueadad|cproductos|ctrx_quarter|ccheques|cinversion|cseguro|tcuentas|vmr_m)"]

print(cols_monetarios)
if (PARAM$RatiosEpico$run) {
  print("procesando RatiosEpico")
  OUTPUT$RatiosEpico$ncol_antes <- ncol(dataset)
  RatiosEpico(dataset, cols_monetarios)
  OUTPUT$RatiosEpico$ncol_despues <- ncol(dataset)
  GrabarOutput()
}

if (PARAM$RatiosEpicoDiv0NA$run) {
  print("procesando RatiosEpico (si es 0 lo hace NA)")
  OUTPUT$RatiosEpicoDiv0NA$ncol_antes <- ncol(dataset)
  RatiosEpicoDiv0HaceNA(dataset, cols_monetarios)
  OUTPUT$RatiosEpicoDiv0NA$ncol_despues <- ncol(dataset)
  GrabarOutput()
}

# esta linea es necesaria?? 
cols_lagueables <- intersect(cols_lagueables, colnames(dataset))

if (PARAM$lag1) {
  print("procesando lag1")
  # creo los campos lags de orden 1
  OUTPUT$lag1$ncol_antes <- ncol(dataset)
  dataset[, paste0(cols_lagueables, "_lag1") := shift(.SD, 1, NA, "lag"),
    by = eval( PARAM$dataset_metadata$entity_id),
    .SDcols = cols_lagueables
  ]

  # agrego los delta lags de orden 1
  for (vcol in cols_lagueables)
  {
    dataset[, paste0(vcol, "_delta1") := get(vcol) - get(paste0(vcol, "_lag1"))]
  }

  OUTPUT$lag1$ncol_despues <- ncol(dataset)
  GrabarOutput()
}


cols_lagueables <- intersect(cols_lagueables, colnames(dataset))
if (PARAM$lag2) {
  print("procesando lag2")
  # creo los campos lags de orden 2
  OUTPUT$lag2$ncol_antes <- ncol(dataset)
  dataset[, paste0(cols_lagueables, "_lag2") := shift(.SD, 2, NA, "lag"),
    by = eval(PARAM$dataset_metadata$entity_id),
    .SDcols = cols_lagueables
  ]

  # agrego los delta lags de orden 2
  for (vcol in cols_lagueables)
  {
    dataset[, paste0(vcol, "_delta2") := get(vcol) - get(paste0(vcol, "_lag2"))]
  }

  OUTPUT$lag2$ncol_despues <- ncol(dataset)
  GrabarOutput()
}


cols_lagueables <- intersect(cols_lagueables, colnames(dataset))
if (PARAM$lag3) {
  print("procesando lag3")
  # creo los campos lags de orden 3
  OUTPUT$lag3$ncol_antes <- ncol(dataset)
  dataset[, paste0(cols_lagueables, "_lag3") := shift(.SD, 3, NA, "lag"),
    by = eval(PARAM$dataset_metadata$entity_id),
    .SDcols = cols_lagueables
  ]

  # agrego los delta lags de orden 3
  for (vcol in cols_lagueables)
  {
    dataset[, paste0(vcol, "_delta3") := get(vcol) - get(paste0(vcol, "_lag3"))]
  }

  OUTPUT$lag3$ncol_despues <- ncol(dataset)
  GrabarOutput()
}

cols_lagueables <- intersect(cols_lagueables, colnames(dataset))
if (PARAM$lag4) {
  print("procesando lag4")
  # creo los campos lags de orden 4
  OUTPUT$lag4$ncol_antes <- ncol(dataset)
  dataset[, paste0(cols_lagueables, "_lag4") := shift(.SD, 4, NA, "lag"),
    by = eval(PARAM$dataset_metadata$entity_id),
    .SDcols = cols_lagueables
  ]

  # agrego los delta lags de orden 4
  for (vcol in cols_lagueables)
  {
    dataset[, paste0(vcol, "_delta4") := get(vcol) - get(paste0(vcol, "_lag4"))]
  }

  OUTPUT$lag4$ncol_despues <- ncol(dataset)
  GrabarOutput()
}


cols_lagueables <- intersect(cols_lagueables, colnames(dataset))
if (PARAM$lag5) {
  print("procesando lag5")
  # creo los campos lags de orden 5
  OUTPUT$lag5$ncol_antes <- ncol(dataset)
  dataset[, paste0(cols_lagueables, "_lag5") := shift(.SD, 5, NA, "lag"),
    by = eval(PARAM$dataset_metadata$entity_id),
    .SDcols = cols_lagueables
  ]

  # agrego los delta lags de orden 5
  for (vcol in cols_lagueables)
  {
    dataset[, paste0(vcol, "_delta5") := get(vcol) - get(paste0(vcol, "_lag5"))]
  }

  OUTPUT$lag5$ncol_despues <- ncol(dataset)
  GrabarOutput()
}



cols_lagueables <- intersect(cols_lagueables, colnames(dataset))
if (PARAM$lag6) {
  print("procesando lag6")
  # creo los campos lags de orden 6
  OUTPUT$lag6$ncol_antes <- ncol(dataset)
  dataset[, paste0(cols_lagueables, "_lag6") := shift(.SD, 6, NA, "lag"),
    by = eval(PARAM$dataset_metadata$entity_id),
    .SDcols = cols_lagueables
  ]

  # agrego los delta lags de orden 6
  for (vcol in cols_lagueables)
  {
    dataset[, paste0(vcol, "_delta6") := get(vcol) - get(paste0(vcol, "_lag6"))]
  }

  OUTPUT$lag6$ncol_despues <- ncol(dataset)
  GrabarOutput()
}



cols_lagueables <- intersect(cols_lagueables, colnames(dataset))
if (PARAM$lag7) {
  print("procesando lag7")
  # creo los campos lags de orden 7
  OUTPUT$lag7$ncol_antes <- ncol(dataset)
  dataset[, paste0(cols_lagueables, "_lag7") := shift(.SD, 7, NA, "lag"),
    by = eval(PARAM$dataset_metadata$entity_id),
    .SDcols = cols_lagueables
  ]

  # agrego los delta lags de orden 7
  for (vcol in cols_lagueables)
  {
    dataset[, paste0(vcol, "_delta7") := get(vcol) - get(paste0(vcol, "_lag7"))]
  }

  OUTPUT$lag7$ncol_despues <- ncol(dataset)
  GrabarOutput()
}

cols_lagueables <- intersect(cols_lagueables, colnames(dataset))
if (PARAM$lag8) {
  print("procesando lag8")
  # creo los campos lags de orden 7
  OUTPUT$lag8$ncol_antes <- ncol(dataset)
  dataset[, paste0(cols_lagueables, "_lag8") := shift(.SD, 8, NA, "lag"),
    by = eval(PARAM$dataset_metadata$entity_id),
    .SDcols = cols_lagueables
  ]

  # agrego los delta lags de orden 8
  for (vcol in cols_lagueables)
  {
    dataset[, paste0(vcol, "_delta8") := get(vcol) - get(paste0(vcol, "_lag8"))]
  }

  OUTPUT$lag8$ncol_despues <- ncol(dataset)
  GrabarOutput()
}

cols_lagueables <- intersect(cols_lagueables, colnames(dataset))
if (PARAM$lag9) {
  print("procesando lag9")
  # creo los campos lags de orden 9
  OUTPUT$lag9$ncol_antes <- ncol(dataset)
  dataset[, paste0(cols_lagueables, "_lag9") := shift(.SD, 9, NA, "lag"),
    by = eval(PARAM$dataset_metadata$entity_id),
    .SDcols = cols_lagueables
  ]

  # agrego los delta lags de orden 9
  for (vcol in cols_lagueables)
  {
    dataset[, paste0(vcol, "_delta9") := get(vcol) - get(paste0(vcol, "_lag9"))]
  }

  OUTPUT$lag9$ncol_despues <- ncol(dataset)
  GrabarOutput()
}


#--------------------------------------
# agrego las tendencias

# ordeno el dataset por <numero_de_cliente, foto_mes> para poder hacer lags
#  es MUY  importante esta linea
setorderv(dataset, PARAM$dataset_metadata$primarykey)

cols_lagueables <- intersect(cols_lagueables, colnames(dataset))




if (PARAM$Tendencias1$run) {
  print("procesando tendencias1")
  OUTPUT$TendenciasYmuchomas1$ncol_antes <- ncol(dataset)
  dataset <- TendenciaYmuchomas(dataset,
    cols = cols_lagueables,
    ventana = PARAM$Tendencias1$ventana, # 6 meses de historia
    tendencia = PARAM$Tendencias1$tendencia,
    minimo = PARAM$Tendencias1$minimo,
    maximo = PARAM$Tendencias1$maximo,
    promedio = PARAM$Tendencias1$promedio,
    ratioavg = PARAM$Tendencias1$ratioavg,
    ratiomax = PARAM$Tendencias1$ratiomax,
    autocorr = PARAM$Tendencias1$autocorr
  )

  OUTPUT$TendenciasYmuchomas1$ncol_despues <- ncol(dataset)
  GrabarOutput()
}


cols_lagueables <- intersect(cols_lagueables, colnames(dataset))
if (PARAM$Tendencias2$run) {
  print("procesando tendencias2")
  OUTPUT$TendenciasYmuchomas2$ncol_antes <- ncol(dataset)
  dataset <- TendenciaYmuchomas(dataset,
    cols = cols_lagueables,
    ventana = PARAM$Tendencias2$ventana, # 6 meses de historia
    tendencia = PARAM$Tendencias2$tendencia,
    minimo = PARAM$Tendencias2$minimo,
    maximo = PARAM$Tendencias2$maximo,
    promedio = PARAM$Tendencias2$promedio,
    ratioavg = PARAM$Tendencias2$ratioavg,
    ratiomax = PARAM$Tendencias2$ratiomax,
    autocorr = PARAM$Tendencias2$autocorr
  )

  OUTPUT$TendenciasYmuchomas2$ncol_despues <- ncol(dataset)
  GrabarOutput()
}

if (PARAM$CanaritosAsesinos1$ratio > 0.0) {
  print("procesando CanaritosAsesinos1 - Pre random forest")
  OUTPUT$CanaritosAsesinos1$ncol_antes <- ncol(dataset)
  CanaritosAsesinos(
    canaritos_ratio = PARAM$CanaritosAsesinos1$ratio,
    canaritos_desvios = PARAM$CanaritosAsesinos1$desvios,
    canaritos_semilla = PARAM$CanaritosAsesinos1$semilla
  )

  OUTPUT$CanaritosAsesinos1$ncol_despues <- ncol(dataset)
  GrabarOutput()
}

#------------------------------------------------------------------------------
# Agrego variables a partir de las hojas de un Random Forest 

if (PARAM$RandomForest$run) {
  print("procesando RandomForest")
  OUTPUT$AgregaVarRandomForest$ncol_antes <- ncol(dataset)
  AgregaVarRandomForest(
    num.trees = PARAM$RandomForest$num.trees,
    max.depth = PARAM$RandomForest$max.depth,
    min.node.size = PARAM$RandomForest$min.node.size,
    mtry = PARAM$RandomForest$mtry,
    semilla = PARAM$RandomForest$semilla
  )

  OUTPUT$AgregaVarRandomForest$ncol_despues <- ncol(dataset)
  GrabarOutput()
  gc()
}


#--------------------------------------------------------------------------
# Elimino las variables que no son tan importantes en el dataset
# with great power comes grest responsability
if (PARAM$CanaritosAsesinos2$ratio > 0.0) {
  print("procesando CanaritosAsesinos2 - final")
  OUTPUT$CanaritosAsesinos2$ncol_antes <- ncol(dataset)
  CanaritosAsesinos(
    canaritos_ratio = PARAM$CanaritosAsesinos2$ratio,
    canaritos_desvios = PARAM$CanaritosAsesinos2$desvios,
    canaritos_semilla = PARAM$CanaritosAsesinos2$semilla
  )

  OUTPUT$CanaritosAsesinos2$ncol_despues <- ncol(dataset)
  GrabarOutput()
}


print("Guardando resultados")
#------------------------------------------------------------------------------
# grabo el dataset

fwrite(dataset,
  file = "dataset.csv.gz",
  logical01 = TRUE,
  sep = ","
)

# copia la metadata sin modificar
write_yaml( PARAM$dataset_metadata, 
  file="dataset_metadata.yml" )

#------------------------------------------------------------------------------

# guardo los campos que tiene el dataset
tb_campos <- as.data.table(list(
  "pos" = 1:ncol(dataset),
  "campo" = names(sapply(dataset, class)),
  "tipo" = sapply(dataset, class),
  "nulos" = sapply(dataset, function(x) {
    sum(is.na(x))
  }),
  "ceros" = sapply(dataset, function(x) {
    sum(x == 0, na.rm = TRUE)
  })
))

fwrite(tb_campos,
  file = "dataset.campos.txt",
  sep = "\t"
)

#------------------------------------------------------------------------------
OUTPUT$dataset$ncol <- ncol(dataset)
OUTPUT$dataset$nrow <- nrow(dataset)

OUTPUT$time$end <- format(Sys.time(), "%Y%m%d %H%M%S")
GrabarOutput()

# dejo la marca final
cat(format(Sys.time(), "%Y%m%d %H%M%S"), "\n",
  file = "z-Rend.txt",
  append = TRUE
)

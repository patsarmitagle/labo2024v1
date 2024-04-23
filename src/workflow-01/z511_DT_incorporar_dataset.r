# Experimentos Colaborativos Default
# Workflow  Catastrophe Analysis

# limpio la memoria
rm(list = ls(all.names = TRUE)) # remove all objects
gc(full = TRUE) # garbage collection

require("data.table")
require("yaml")

# Parametros del script
PARAM <- read_yaml( "parametros.yml" )

# FIN Parametros del script

OUTPUT <- list()

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

GrabarOutput <- function() {
  write_yaml(OUTPUT, file = "output.yml") # grabo output
}
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Aqui empieza el programa
OUTPUT$PARAM <- PARAM
OUTPUT$time$start <- format(Sys.time(), "%Y%m%d %H%M%S")


# cargo el dataset
dataset <- fread( paste0( "./datasets/", PARAM$archivo) )

#--------------------------------------
# verifico que existan los campos de la metadata

campos <- copy( colnames( dataset ) )

campitos <- c( PARAM$primarykey, PARAM$entity_id, PARAM$periodo, PARAM$clase )
campitos <- unique( campitos )

for( vcampo in campitos )
  if( ! (vcampo %in% campos ) ) 
  {
    t <- format(Sys.time(), "%Y%m%d %H%M%S")
    cat( t, "\n",
      file = "z-Rabort.txt",
      append = TRUE
    )

    cat( t, "\n",
      file = "z-Rabort-hist.txt",
      append = TRUE
    )

    cat( "No existe el campo : ", vcampo, "\n" )

    stop( "Error faltal" )
  }

GrabarOutput()
#--------------------------------------
# verifico primarykey

pk_qty <- nrow(unique(dataset[ ,PARAM$primarykey, with=FALSE ]))

if( pk_qty != nrow( dataset ) )
{
  t <- format(Sys.time(), "%Y%m%d %H%M%S")
  cat( t, "\n",
    file = "z-Rabort.txt",
    append = TRUE
  )

  cat( t, "\n",
    file = "z-Rabort-hist.txt",
    append = TRUE
  )

  cat( "Primary Key inconsistente" )

  stop( "Error faltal" )
  }

#--------------------------------------

dataset_metadata <- copy( PARAM )
dataset_metadata$archivo <- NULL
dataset_metadata$semilla <- NULL

write_yaml( dataset_metadata, file="dataset_metadata.yml" )

#--------------------------------------
# ordeno el dataset

setorderv(dataset, PARAM$primarykey)

#------------------------------------------------------------------------------
# grabo el dataset

fwrite(dataset,
  file = "dataset.csv.gz",
  logical01 = TRUE,
  sep = ","
)

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

#------------------------------------------------------------------------------
# dejo la marca final

cat(format(Sys.time(), "%Y%m%d %H%M%S"), "\n",
  file = "z-Rend.txt",
  append = TRUE
)

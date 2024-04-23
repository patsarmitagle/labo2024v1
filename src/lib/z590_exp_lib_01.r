# libreria complicada, que accede a variables de sus padres
#   para simplificar la vida a quien la usa

#------------------------------------------------------------------------------
# corre en Linux & friends

exp_softlink <- function( dest, source)
{
  st <- paste0( "ln -sf ", source, "   ", dest )
  system( st )
}
#------------------------------------------------------------------------------

exp_wf_init <- function( pnombrewf, pvirgen=FALSE )
{
  envg$EXPENV$nombrewf <- pnombrewf
  
  if( !exists( "pnombrewf" ) )
    stop( "Falla catastrofica en exp_wf_init() todo workflow debe tener un nombre." )

  setwd( envg$EXPENV$wf_dir )
  existe <- dir.exists( pnombrewf )

  if( pvirgen==TRUE & existe==TRUE ) {
    cat( format(Sys.time(), "%Y%m%d %H%M%S"),
      pnombrewf,
      "\n",
      sep="\t",
      file = "z-Rabort.txt",
      append = TRUE 
    )

    stop( "Falla catastrofica, ya existe un workflow con el mismo nombre  y el llamado fue virgen=TRUE exp_wf_init() \n" )
  }

  dir.create( pnombrewf, showWarnings = FALSE)
  setwd( pnombrewf )

  resultado <- 0
  # if( file.exists("z-Rend.txt") )  return(-1L)

  return( resultado )
}
#------------------------------------------------------------------------------

exp_wf_end <- function( pnombrewf, pvirgen=FALSE )
{
  if( !exists( "pnombrewf" ) )
    stop( "Falla catastrofica, todo workflow debe tener un nombre." )

  setwd( envg$EXPENV$wf_dir )

  # no puede suceder esto en este punto
  if( !dir.exists( pnombrewf ) ) {
    stop( paste0( "Falla catastrofica, al final workflow pero no existe la carpeta", pnombrewf) )
  }

  # agrego que termino todo OK
  setwd( pnombrewf )
  # dejo la marca final
  cat( format(Sys.time(), "%Y%m%d %H%M%S"), "\n",
    file = "z-Rend.txt",
    append = TRUE
  )

  # vuelvo a la carpeta local de workflows
  setwd( envg$EXPENV$wf_dir_local )
}
#------------------------------------------------------------------------------

exp_init_datos <- function( pmyexp, parch, pserver )
{
  if( !exists( "pmyexp" ) )
    stop( "Falla catastrofica  en  exp_init_datos()  todo experimento debe tener un nombre." )

  if( !exists( "parch" ) )
    stop( "Falla catastrofica  en  exp_init_datos()  debe haber nombre de archivo." )

  if( !file.exists( paste0(envg$EXPENV$datasets_dir, parch) ) )
    stop( paste0("Falla catastrofica, no existe el archivo ", parch) )

  output <- list()
  output$resultado <- 0L
  output$meta$myexp <- pmyexp
  output$meta$server <- pserver
  output$archivo <- parch
  output$gcloud <- envg$EXPENV$gcloud

  setwd( envg$EXPENV$exp_dir )
  dir.create( pmyexp, showWarnings = FALSE)
  setwd( pmyexp )

  # si ya esta todo procesado de antes
  if( file.exists("z-Rend.txt") ) {
    output$resultado <- -1L
    return(output)
  }


  # horrible, traigo el nombre del workflow, variable el environment previo
  env_up <- parent.frame( n=2 )
  nombrewf <- env_up$pnombrewf
  cat( "nombrewf = ", nombrewf, "\n" )
  setwd( envg$EXPENV$wf_dir )
  setwd( nombrewf )
  if( !dir.exists(pmyexp) )
    exp_softlink( pmyexp, paste0( envg$EXPENV$exp_dir, "/", pmyexp) )

  # unico link en este caso, a la carpeta  datasets
  setwd( envg$EXPENV$exp_dir )
  setwd( pmyexp )
  if( !dir.exists("datasets") )
    exp_softlink( "datasets", envg$EXPENV$datasets_dir )

  # quede en la carpeta del experimento
  return(output)
}
#------------------------------------------------------------------------------

exp_init <- function( pmyexp, pinputexps, pserver )
{
  if( !exists( "pmyexp" ) )
    stop( "Falla catastrofica, todo experimento debe tener un nombre." )

  if( !exists( "pinputexps" ) )
    stop( "Falla catastrofica, debe existir un vector de experimentos input." )

 
  output <- list()
  output$resultado <- 0L
  output$meta$myexp <- pmyexp
  output$meta$server <- pserver
  output$input <- pinputexps

  output$gcloud <- envg$EXPENV$gcloud

  setwd( envg$EXPENV$exp_dir )
  dir.create( pmyexp, showWarnings = FALSE)
  setwd( pmyexp )

  # si ya esta todo procesado de antes
  if( file.exists("z-Rend.txt") ) {
    output$resultado <- -1L
    return(output)
  }

  setwd( envg$EXPENV$exp_dir )
  for( vexp in pinputexps )
  {
    if( !dir.exists( vexp ) )
    {
      cat( "Fatal Error al querer correr ",  pmyexp, "\n" )
      cat( "NO existe el experimento input :",  vexp, "\n" )
      stop("Saliendo")
    }
  }

  # creo los softlinks internos a experimentos input
  setwd( envg$EXPENV$exp_dir )
  setwd( pmyexp )
  for( vexp in pinputexps )
  {
    if( !dir.exists(vexp) )
      exp_softlink( vexp, paste0( envg$EXPENV$exp_dir, "/", vexp) )
  }

  # creo el link al exp dentro del workflow
  # horrible, traigo el nombre del workflow, variable el environment previo
  env_up <- parent.frame( n=2 )
  nombrewf <- env_up$pnombrewf

  setwd( envg$EXPENV$wf_dir )
  setwd( nombrewf )
  if( !dir.exists(pmyexp) )
    exp_softlink( pmyexp, paste0( envg$EXPENV$exp_dir, "/", pmyexp) )

  # quedo en la carpeta del experimento
  setwd( envg$EXPENV$exp_dir )
  setwd( pmyexp )

  return(output)
}
#------------------------------------------------------------------------------

exp_end <- function( pmyexp )
{
  if( !exists( "pmyexp" ) )
    stop( "Falla catastrofica en exp_end_() todo experimento debe tener un nombre." )

  env_up <- parent.frame( n=3 )
  nombrewf <- env_up$pnombrewf
  cat( "exp_end  nombrewf:", nombrewf, "\n" )

  setwd( envg$EXPENV$exp_dir )
  if( !dir.exists( pmyexp ) )
    stop( paste0( "Falla catastrofica, al final el experimento no existe la carpeta exp", pmyexp ))

  setwd( envg$EXPENV$wf_dir )
  setwd( nombrewf )
  if( !dir.exists( pmyexp ) )
    stop( paste0( "Falla catastrofica, al final el experimento no existe la carpeta en el workflow", pmyexp ))

  setwd( envg$EXPENV$exp_dir )
  dir.create( pmyexp, showWarnings = FALSE)
  setwd( pmyexp )

  # dejo la marca final
  cat(format(Sys.time(), "%Y%m%d %H%M%S"), "\n",
  file = "z-Rend.txt",
  append = TRUE
  )

}
#------------------------------------------------------------------------------

exp_extraer_archivo <- function( st )
{
  largo <- nchar(st)
  i <- largo
  while( i >= 1 & substr( st, i, i) != "/" )  i <- i - 1

  if( substr(st, i, i) == "/" ) i <- i + 1

  res <- ""
  if( i <= largo )  res <- substr( st, i, largo )

  return( res )
}
#------------------------------------------------------------------------------

exp_correr_script <- function( pparam_local )
{
  if( !exists( "pparam_local" ) )
    stop( "Falla catastrofica, se le debe pasar parametro a  exp_correr_script()" )

  script_origen <- paste0( envg$EXPENV$repo_dir, "/", pparam_local$meta$script )
  if( !file.exists( script_origen ) )
    stop( paste0("Falla catastrofica, no exists el script", script_origen) )

  # se supone que estamos en la carepta
  # grabo el archivo de parametros que va a recibir el script local
  param_simple <- data.table::copy( pparam_local )
  param_simple$experimento <- param_simple$meta$myexp
  param_simple$resultado <- NULL
  param_simple$meta <- NULL
  param_simple$gcloud <- NULL

  arch_semillas <- paste0( envg$EXPENV$datasets_dir, envg$EXPENV$arch_sem)
  tabla_semillas <- fread( arch_semillas )
  ksemilla_azar <- tabla_semillas[ 1, semilla ]  # 1 es mi primera semilla
  param_simple$semilla <- ksemilla_azar
  write_yaml( param_simple, "parametros.yml" )


  script_corto <- exp_extraer_archivo( pparam_local$meta$script )
  # copio el script R a la carpeta del experimento
  #   esto deberia manejarse con Git
  if( !file.exists(script_corto)) {
    file.copy( script_origen, "." )
    Sys.chmod( script_corto, mode = "544", use_umask = TRUE)
  }

  # creo el script que corre el experimento
  if( !file.exists("run.sh"))
  {
    linea1  <- "tabulador=\"\t\"\n"
    linea2  <- "echo \"timestamp\tevent\"  >  log.txt \n"
    linea3  <- "fecha0=$(date +\"%Y%m%d %H%M%S\") \n"
    linea4  <- "echo \"$fecha0\"\"$tabulador\"\"SH_START\" >> log.txt \n"

    linea5 <- paste0( "Rscript --vanilla ",
                      script_corto,
                      "  " ,
                      "parametros.yml",
                      "  2>&1 | tee outfile \n" )

    linea6  <- "fecha1=$(date +\"%Y%m%d %H%M%S\") \n"
    linea7  <- "echo \"$fecha1\"\"$tabulador\"\"SH_END\" >> log.txt \n"

    comando  <- paste0( linea1, linea2, linea3, linea4, linea5, linea6, linea7)

    # creacion del archivo
    shell_script <- "run.sh"
    cat( comando, 
         file= shell_script )

    #doy permisos de ejecucion al shell script
    Sys.chmod( shell_script, mode = "744", use_umask = TRUE)
  }

  if( file.exists( "z-Rabort.txt" ) )  file.remove( "z-Rabort.txt" )

  #ejecuto en linux el script recien creado
  if( pparam_local$meta$server == "local" )
     system( paste0( "./", shell_script ) )

  if( pparam_local$meta$server == "gcloud" )
  {
     # aqui va todo el circo de crear una virtual machine
  }

  # el proceso hijo aborto,  yo tambien debo abortar
  if( file.exists( "z-Rabort.txt" ) )
    stop( paste0("Falla catastrofica, ha terminado anormalmente", script_origen) )

  exp_end( pparam_local$meta$myexp )
  return( pparam_local$meta$myexp  )
}
#------------------------------------------------------------------------------

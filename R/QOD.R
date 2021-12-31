#' A quality of data inspector
#'
#' @description   a QoD inspector class
#' @export
QOD <- function( UM = "" ) {

  global.dataLoader <- c();
  global.UM <- ""
  
  #===========================================================
  # loadDataset
  #===========================================================
  loadDataset<-function( dataList ) {
    # Clear some possible previously inserted attribute
    clearAttributes()
    # set the new attributes
    EFOMM <- EfirstOrderMarkovModel()
    EFOMM$loadDataset( dataList )
    global.EFPTs <<- EFOMM$getEFPTs()
    # set the new attributes
    global.dataLoader <<- dataList
  }
  #=================================================================================
  # query
  #=================================================================================
  query <- function( from , to , complement = FALSE, time.range=c(0,Inf), step.range = c(1,Inf) , UM = NA, 
                     arr.passingThrough = c(), arr.NOTpassingThrough = c(),
                     forceCheck = TRUE) {
    EventName <- global.dataLoader$csv.EVENTName
    if( (!(from %in% global.dataLoader$arrayAssociativo) | !(to %in% global.dataLoader$arrayAssociativo)) & forceCheck == TRUE) {
      stop("Error: from or to not available as events in the Event Log")
    }
    
    if(is.na(UM)) UM <- global.UM
    tmp.res <- unlist(lapply( names( global.dataLoader$pat.process ), function(ID) {
      subMM <- global.dataLoader$pat.process[[ID]]
      
      if( UM == "hours" ) subMM$pMineR.deltaDate <- subMM$pMineR.deltaDate/60
      if( UM == "days" ) subMM$pMineR.deltaDate <- subMM$pMineR.deltaDate/(60*24)
      if( UM == "weeks" ) subMM$pMineR.deltaDate <- subMM$pMineR.deltaDate/(60*24*7)
      if( UM == "months" ) subMM$pMineR.deltaDate <- subMM$pMineR.deltaDate/(60*24*7*30)
      
      begin.line <- subMM[1,]; begin.line[ EventName ] <- "BEGIN"
      end.line <- subMM[nrow(subMM),]; end.line[ EventName ] <- "END"
      subMM <- rbind( begin.line, subMM ); subMM <- rbind( subMM, end.line )

      arr.from <- which(subMM[[EventName]] == from)
      arr.to <- which(subMM[[EventName]] == to)
      
      if( length(arr.from) == 0 | length(arr.to) == 0 ) return( FALSE )
      
      MM <- expand.grid.unique(arr.from,arr.to)
      MM <- cbind(MM , rep(0,nrow(MM))); MM <- cbind(MM , rep(0,nrow(MM))); MM <- cbind(MM , rep(0,nrow(MM)))
      tmp <- lapply(1:nrow(MM),function(riga) { 
        MM[riga,3] <<- MM[riga,2] - MM[riga,1]
        MM[riga,4] <<- subMM$pMineR.deltaDate[ MM[riga,2] ] - subMM$pMineR.deltaDate[ MM[riga,1] ]
      })
      colnames(MM) <- c("from","to","step","time","valid")
      
      for( riga in 1:nrow(MM)) {
        valido <- TRUE
        if(! (MM[riga,"step"] >= step.range[1] & MM[riga,"step"] <= step.range[2]) ) {
          valido <- FALSE
        }
        if(! (MM[riga,"time"] >= time.range[1] & MM[riga,"time"] <= time.range[2]) ) {
          valido <- FALSE
        }
        if( sum(arr.passingThrough %in% subMM[[EventName]][ MM[riga,1]:MM[riga,2] ]) != length(arr.passingThrough) )  {
          valido <- FALSE
        }
        if( sum(arr.NOTpassingThrough %in% subMM[[EventName]][ MM[riga,1]:MM[riga,2] ]) >0 )  {
          valido <- FALSE
        }        
        # browser()
        MM[riga,"valid"] <- valido
      }
      return(sum(MM[,"valid"]) > 0)
    } ))
    
    res <- names( global.dataLoader$pat.process )[which(tmp.res==TRUE)]
    
    if( complement == TRUE ) { 
      res <- names(global.dataLoader$pat.process)[which( !(names(global.dataLoader$pat.process) %in% res))]
    }
    
    return(res)
  } 
  
  #=================================================================================
  # clearAttributes
  #=================================================================================
  clearAttributes<-function() {
    costructor( UM = UM );
  }  
  #===========================================================
  # costructor
  # E' il costruttore della classe
  #===========================================================
  costructor<-function( UM ) {
    global.dataLoader <<- ''
    global.EFPTs <<- c()
    global.UM <<- UM
    if( global.UM == "" ) global.UM <<- "days"
  }
  #===========================================================
  costructor( UM = UM );
  #===========================================================
  return( list(
    "loadDataset"=loadDataset,
    "query"=query
  ) )
}

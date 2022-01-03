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
                     forceCheck = TRUE, returnCompleteMatrix = FALSE) {
    EventName <- global.dataLoader$csv.EVENTName
    if( (!(from %in% global.dataLoader$arrayAssociativo) | !(to %in% global.dataLoader$arrayAssociativo)) & forceCheck == TRUE) {
      stop("Error: from or to not available as events in the Event Log")
    }
    
    if(is.na(UM)) UM <- global.UM
    mainMM <- c()
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
      mainMM <<- rbind(mainMM, cbind( rep(ID,nrow(MM)), MM ) )
      return(sum(MM[,"valid"]) > 0)
    } ))
    # browser()
    # mainMM <- mainMM[ which(mainMM[,"valid"] == 1),]
    nomiColonne <- colnames(mainMM)
    mainMM <- matrix(mainMM[ which(mainMM[,"valid"] == 1),], ncol=ncol(mainMM))
    if(length(mainMM) == 0) return(NA)
    colnames(mainMM) <- nomiColonne
    mainMM[,2] <- as.numeric(mainMM[,2])-1
    mainMM[,3] <- as.numeric(mainMM[,3])-1
    if( returnCompleteMatrix == TRUE ) return(mainMM)
    
    res <- names( global.dataLoader$pat.process )[which(tmp.res==TRUE)]
    
    if( complement == TRUE ) { 
      res <- names(global.dataLoader$pat.process)[which( !(names(global.dataLoader$pat.process) %in% res))]
    }
    
    return(res)
  } 
  #=================================================================================
  # path.count
  #=================================================================================  
  path.count <- function( evtSequenceLength = 1, notMoreThanOnePerPatient = FALSE , fromBegin = TRUE ) {
    browser()
  }
  #=================================================================================
  # eventHeatmap
  #=================================================================================
  eventHeatmap <- function( cex = 0.5 , threshold.low = 0.5, threshold.hi = 1, show.diagonal = TRUE, par.margin = c(4, 10, 10, 2)) {
    
    objDL.new.export <- global.dataLoader
    
    arr.eventi <- objDL.new.export$arrayAssociativo[!(objDL.new.export$arrayAssociativo %in% c("BEGIN","END"))]
    MM.Cross <- matrix( 0,nrow = length(arr.eventi), ncol = length(arr.eventi) )
    colnames(MM.Cross) <- arr.eventi; rownames(MM.Cross) <- arr.eventi;
    tmp.1 <- lapply( rownames(MM.Cross) , function(event.C) {
      tmp.2 <- lapply(colnames(MM.Cross), function(event.R) {
        tmp.3 <- lapply( names(objDL.new.export$pat.process) , function(patID) {
          arr.evt.to.chech <- objDL.new.export$pat.process[[patID]][[objDL.new.export$csv.EVENTName]]
          if( event.C %in% arr.evt.to.chech & event.R %in% arr.evt.to.chech) {
            MM.Cross[ event.R , event.C ] <<- MM.Cross[ event.R , event.C ] + 1
          }
        })
      } )
    })
    aaa <- MM.Cross
    tmp.1 <- lapply( 1:nrow(MM.Cross) , function(riga) { 
      MM.Cross[riga,] <<- MM.Cross[riga,] / MM.Cross[riga,riga]
    })
    
    
    par(mar = par.margin) 
    image(t(MM.Cross[nrow(MM.Cross):1,]),col = heat.colors(256) , axes=FALSE )
    arr.posizioni <- (0.1:ncol(MM.Cross)/(ncol(MM.Cross)-1))
    axis(2, arr.posizioni, labels=rownames(MM.Cross)[length(rownames(MM.Cross)):1],las=2)
    axis(3, arr.posizioni, labels=rownames(MM.Cross),las=2)
    for( riga in 1:nrow(MM.Cross)) {
      for( colonna in 1:ncol(MM.Cross)) {
        valore <- t(MM.Cross[ncol(MM.Cross)-colonna+1,riga])
        if( valore >= threshold.low & valore <= threshold.hi ) {
          text((riga-1)/(nrow(MM.Cross)-1),(colonna-1)/(ncol(MM.Cross)-1),format(valore,digits = 2) , cex=cex ) 
        }
      }
    }
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
    "query"=query,
    "eventHeatmap"=eventHeatmap
    # "path.count"=path.count
  ) )
}

# =============================================================================
#' Some useful tools
#' 
#' @description  A class which provide some tools. pMineR intarnal use only.
# =============================================================================
utils<-function() {
  dectobin <- function(y) {
    # find the binary sequence corresponding to the decimal number 'y'
    stopifnot(length(y) == 1, mode(y) == 'numeric')
    q1 <- (y / 2) %/% 1
    r <- y - q1 * 2
    res = c(r)
    while (q1 >= 1) {
      q2 <- (q1 / 2) %/% 1
      r <- q1 - q2 * 2
      q1 <- q2
      res = c(r, res)
    }
    return(res)
  } 
  is.included<-function( a , b ) {
    if(sum(is.element(a,b)) == length(a)) return(TRUE)
    else return(FALSE)
  }
  cleanUTF <- function( dati , colonna.evento , def.val.to.substitute = 95 ){
    for (riga in 1:nrow(dati)){
      arr <- as.numeric(charToRaw(x = as.character(dati[riga,colonna.evento]) ))
      arr[ which(arr > 127)]<-def.val.to.substitute 
      dati[riga,colonna.evento] <- intToUtf8(arr)
    }
    return(dati)
  }  
  format.data.for.csv<-function(listaProcessi, lista.validi, typeOfRandomDataGenerator="dayAfterDay",output.format.date = "%d/%m/%Y %H:%M:%S") { 
    big.csv<-c()
    ct <- 1
    
    for(i in names(listaProcessi)) {
      numeroElementi<-length(listaProcessi[[i]])
      
      if(typeOfRandomDataGenerator=="dayAfterDay") giorni.da.sommare <- as.integer(runif(n = numeroElementi,min=1,max=1))
      if(typeOfRandomDataGenerator=="randomDay1-4") giorni.da.sommare <- as.integer(runif(n = numeroElementi,min=1,max=4) )
      if(typeOfRandomDataGenerator=="randomWeek1-4") giorni.da.sommare <- as.integer(runif(n = numeroElementi,min=1,max=4) * 7)
      if(typeOfRandomDataGenerator=="randomMonth1-4") giorni.da.sommare <- as.integer(runif(n = numeroElementi,min=1,max=4) * 30)
      
      array.Date <- as.character(format(as.Date("01/01/2000 12:00:00",format=output.format.date) + cumsum(giorni.da.sommare) ,format=output.format.date) )
      matrice<-cbind(rep(ct,numeroElementi),listaProcessi[[i]],array.Date,rep(as.character(lista.validi[ct]),numeroElementi) )
      big.csv<-rbind(big.csv,matrice )
      ct <- ct + 1
    }
    if(!is.null(dim(big.csv))) {
      colnames(big.csv)<-c("patID","event","date","valido")
    }
    return(big.csv)
  }  
  # check if the file is pure ASCII
  # 
  # a simple function able to check if the indicated file is a pure ASCII file
  # the name of the file that need to be checked  
  IsASCII<-function( fileName ) {
    return(c_IsASCII(fileName = fileName))
  }  
  return(list(
    "dectobin" = dectobin,
    "is.included" = is.included,
    "format.data.for.csv" = format.data.for.csv,
    "cleanUTF"=cleanUTF,
    "IsASCII"=IsASCII
  ))
}
# =============================================================================
#' textObj
#' Una classe ad uso interno per manipolare testi
# =============================================================================
textObj<-function() {
  testo<-'';
  add<-function( stringa, carriage=TRUE) {
    if(length(stringa)>1) stringa<-paste(stringa,collapse='')
    if(carriage==TRUE)
      testo <<- paste( c(testo,'\n',stringa), collapse = ''  ) 
    else
      testo <<- paste( c(testo,stringa), collapse = ''  ) 
  }
  get<-function() {
    return(testo)
  } 
  costructor<-function() {
    testo<<-'';
  }
  return(list("add"=add,"get"=get))
}

#' some data processing useful tools
#' 
#' @description  A class which provide some tools. pMineR intarnal use only.
dataProcessor<-function() {
  #=================================================================================
  # buildMMMatrices.and.other.structures
  # costruisce la MM matrix ed anche altra robaccia
  #=================================================================================    
  buildMMMatrices.and.other.structures<-function(mydata, EVENT.list.names, 
                                                 EVENTName, EVENTDateColumnName=NA, 
                                                 ID.act.group,
                                                 max.char.length.label = 50,
                                                 verbose.mode = TRUE) {
    
    # costruisci la matrice
    MM<-matrix(0, ncol=length(unique(mydata[[EVENT.list.names]]))+2, nrow=length(unique(mydata[[EVENT.list.names]]))+2 )
    colnames(MM)<-c("BEGIN","END",unique(as.character(mydata[[EVENT.list.names]])))
    rownames(MM)<-colnames(MM)
    # browser()
    if(("" %in% trimws(colnames(MM))) == TRUE) {
      return( list("error"=TRUE, "errCode"=1)  )
    }
    
    if(max(nchar(colnames(MM)))>max.char.length.label)  {
      return( list("error"=TRUE, "errCode"=2)  )
    }
    if(length(grep("'", colnames(MM))))  {
      return( list("error"=TRUE, "errCode"=3)  )
    }    
    
    # Creiamo anche la matrice con le density dei tempi di transizione
    # (ma solo se c'e' un campo DATA TIME)
    MM.den.list<-list()
    MM.den.list.high.det<-list()
    lst.aaa <- list()
    
    # ora scorri la storia dei singoli pazienti per estrarre le ricorrenze
    # per ogni paziente
    if( verbose.mode == TRUE ) pb <- txtProgressBar(min = 0, max = length(ID.act.group), style = 3)
    for(patID in seq(1,length(ID.act.group))) {
      if( verbose.mode == TRUE ) setTxtProgressBar(pb, patID)
      # su ogni elemento del percorso clinico
      # t e' il "tempo" in senso di "step"
      for(t in seq(1,nrow(ID.act.group[[patID]]))) {
        # vedi se devi legare il BEGIN
        if( t == 1) {
          valore<-MM[ "BEGIN", ID.act.group[[patID]][ t ,EVENT.list.names] ]
          MM[ "BEGIN", ID.act.group[[patID]][ t ,EVENT.list.names] ]<-valore+1
        }
        # vedi se devi legare l'END   
        if( t == nrow(ID.act.group[[patID]])) {
          nomeCampo<-ID.act.group[[patID]][t,EVENT.list.names]
          MM[nomeCampo,"END"]<-MM[nomeCampo,"END"]+1
        }

        # tutti gli altri
        if( t < nrow(ID.act.group[[patID]])) {
          nomeCampo.pre<-ID.act.group[[patID]][t,EVENT.list.names]
          nomeCampo.post<-ID.act.group[[patID]][t+1,EVENT.list.names]
          MM[ nomeCampo.pre, nomeCampo.post ]<-MM[ nomeCampo.pre, nomeCampo.post ]+1
          if(EVENTDateColumnName!='' & ! is.na(EVENTDateColumnName)){
            delta.date<-as.numeric(difftime(as.POSIXct(ID.act.group[[patID]][t+1,EVENTDateColumnName], format = "%d/%m/%Y %H:%M:%S"),as.POSIXct(ID.act.group[[patID]][t,EVENTDateColumnName], format = "%d/%m/%Y %H:%M:%S"),units = 'mins'))
            if(length(MM.den.list[[ nomeCampo.pre]])==0) MM.den.list[[ nomeCampo.pre]]<-list()
            if(length(MM.den.list[[ nomeCampo.pre]][[ nomeCampo.post ]])==0) MM.den.list[[ nomeCampo.pre]][[ nomeCampo.post ]]<-c()
            MM.den.list[[ nomeCampo.pre]][[ nomeCampo.post ]]<-c(MM.den.list[[ nomeCampo.pre]][[ nomeCampo.post ]],delta.date)
          }
        }    
      }
      # invoca il programma in C per estrarre i tempi reciproci fra TUTTI
      iii <- unlist(lapply(ID.act.group[[patID]][,EVENT.list.names] , function(x) which(colnames(MM)==x) ))
      massimo <-max(iii)
      out.MM<-rep( 0 , (massimo)*(massimo) )
      out.delta<-c()
      nuovoOut <- c()
      
      # aaa <- transitionsTime( iii , ID.act.group[[patID]][,"pMineR.deltaDate"], max(iii) );
      
      # -im 
      kkk <- c()
      righe.massime <- nrow(ID.act.group[[patID]])

      tmp.lll <- lapply(1:(righe.massime-1), function(riga){
        kkk <- ID.act.group[[patID]][(riga+1:righe.massime),c(EVENTName,"pMineR.deltaDate")]
        kkk[,c("pMineR.deltaDate")] <- kkk[,c("pMineR.deltaDate")] -   ID.act.group[[patID]][riga,"pMineR.deltaDate"]
        tmp.z <- lapply(unique(kkk[,1]),function(evt.tmp) {
          tempo <- kkk$pMineR.deltaDate[which(kkk[[EVENTName]]==evt.tmp)[1]]
          return(c(evt.tmp,tempo))
        })
        tmp.z <- matrix(unlist(tmp.z),ncol=2,byrow = T )
        tmp.z <- cbind( rep(ID.act.group[[patID]][riga,EVENTName],nrow(tmp.z)), tmp.z )
        return(tmp.z)
      })
      tmp.ll <- do.call(rbind,tmp.lll)
      tmp.ll <- tmp.ll[which(!is.na(tmp.ll[,2])),]

      if( length(tmp.ll)>3) {
        tmptmp <- apply(tmp.ll,MARGIN = 1,function(riga) {
          # browser()
          if(length(MM.den.list.high.det[[ riga[1] ]])==0) MM.den.list.high.det[[ riga[1] ]] <<-list()
          if(length(MM.den.list.high.det[[ riga[1] ]][[ riga[2] ]])==0) MM.den.list.high.det[[ riga[1] ]][[ riga[2] ]] <<- list()
          MM.den.list.high.det[[ riga[1] ]][[ riga[2] ]] <<- c(MM.den.list.high.det[[ riga[1] ]][[ riga[2] ]],as.numeric(riga[3]))
          
        })
        
      }
      # -fm

      # browser()
      # 
      # mm.in <- matrix(c(iii,ID.act.group[[patID]][,"pMineR.deltaDate"]),nrow=2,byrow = T)
      # mm.out <- t(matrix(c(aaa$from,aaa$to,aaa$time),nrow=3,byrow = T))
      # for( riga in seq(1,nrow(mm.out))) {
      #   int.from <-colnames(MM)[mm.out[riga,1]];
      #   int.to <-colnames(MM)[mm.out[riga,2]];
      #   delta.tempo <-mm.out[riga,3];
      #   if(length(MM.den.list.high.det[[ int.from ]])==0) MM.den.list.high.det[[ int.from]]<-list()
      #   if(length(MM.den.list.high.det[[ int.from]][[ int.to ]])==0) MM.den.list.high.det[[ int.from]][[ int.to ]]<-c()
      #   MM.den.list.high.det[[ int.from]][[ int.to ]]<-c(MM.den.list.high.det[[ int.from]][[ int.to ]],delta.tempo)
      # }
    }
    
    for(primo in names(MM.den.list.high.det)) {
      for(secondo in names(MM.den.list.high.det[[primo]])) {
        MM.den.list.high.det[[primo]][[secondo]] <- unlist(MM.den.list.high.det[[primo]][[secondo]])
      }
    }
    # browser()
    if( verbose.mode == TRUE ) close(pb)
    quanti.da.fare<-length(names(MM.den.list)) * length(names(MM.den.list))
    # save(lst.aaa,file = "c://projects/lst.aaa.RData")
    # Calcola la matrice delle medie dei tempi
    # Sarebbe bello avere le density... vabbe'. piu' avanti
    if(EVENTDateColumnName!='' & !is.na(EVENTDateColumnName)){
      MM.mean.time<-MM
      MM.mean.time[ 1:nrow(MM.mean.time) , 1:ncol(MM.mean.time)   ]<-Inf
      for(state.from in names(MM.den.list))  {
        for(state.to in names(MM.den.list[[state.from]]))  {
          MM.mean.time[state.from,state.to ]<-mean(MM.den.list[[ state.from]][[ state.to ]])
        }        
      }
    }
    
    # CALCOLO LA MATRICE DEI FLUSSI FUORI DALLO STATO
    
    if(EVENTDateColumnName!='' & !is.na(EVENTDateColumnName)){
      MM.mean.outflow.time<-MM
      MM.mean.outflow.time[ 1:nrow(MM.mean.outflow.time) , 1:ncol(MM.mean.outflow.time)   ]<-NA
      for(state.from in names(MM.den.list))  {
        for(state.to in names(MM.den.list[[state.from]]))  {
          MM.mean.outflow.time[state.from,state.to ]<-mean(MM.den.list[[ state.from]][[ state.to ]][which(MM.den.list[[ state.from]][[ state.to ]] >=0 & state.from != state.to)])
        }
      }
    }
    
    # costruisci una semplice versione, con le parole (come piace tanto a Van der Aalst)
    wordSequence.TMP01<-list();
    for(i in seq(1,length(ID.act.group))) {
      IDPat<-names(  ID.act.group)[i]
      wordSequence.TMP01[[IDPat]]<-ID.act.group[[ IDPat ]][[EVENTName]]
    }    
    return(list( "arrayAssociativo" = rownames(MM),
                 "footPrint"="",
                 "MMatrix"=MM,
                 "MM.mean.time"=MM.mean.time,
                 "MM.density.list"=MM.den.list,
                 "MM.mean.outflow.time"=MM.mean.outflow.time,
                 "MM.den.list.high.det" = MM.den.list.high.det,
                 "pat.process"=ID.act.group,
                 "wordSequence.raw"=wordSequence.TMP01,
                 "error"=FALSE) )    
  }  
  
  #=================================================================================
  # createSequenceMatrix
  # crea una matrice di transizione a partire da una mera sequenza di eventi.
  # Creata per poter evitare di dover usare il pacchetto markovChain
  #=================================================================================      
  createSequenceMatrix<-function( sequence2parse ) {
    
    sequenza.simboli <- unique(as.character(sequence2parse))
    MM<-matrix(0, ncol=length(sequenza.simboli), nrow=length(sequenza.simboli) )  
    colnames(MM)<-sequenza.simboli
    rownames(MM)<-sequenza.simboli
    
    # cicla su ogni elemento della sequenza ed incrementa la relativa posizione nella 
    # matrice di transizione from=>to
    for(t in seq(1,length(sequence2parse)-1)) {
      # tutti gli altri
      nomeCampo.pre<-sequence2parse[t]
      nomeCampo.post<-sequence2parse[t+1]
      MM[ nomeCampo.pre, nomeCampo.post ]<-MM[ nomeCampo.pre, nomeCampo.post ]+1
    }
    return(list(
      "transitionCountMatrix" = MM
    ))
  }
  
  return(list(
    "buildMMMatrices.and.other.structures"=buildMMMatrices.and.other.structures,
    "createSequenceMatrix" = createSequenceMatrix
  ))
}

#' calcola EFT
#' 
#' @description  Funzione per calcolare la EFM
#' @export
calcolaEnhancedFootPrintTable.pat.process <- function( dataLoaderOBJ , skip.less.than = 0, threshold.perc = 0.00000001 ) {

  EFPT.neverAfter<-dataLoaderOBJ$MMatrix;  EFPT.neverAfter[,]<-"!>>"
  EFPT.alwaysAfter<-dataLoaderOBJ$MMatrix;  EFPT.alwaysAfter[,]<-">>"
  # EFPT.neverBefore<-dataLoaderOBJ$MMatrix;  EFPT.neverBefore[,]<-"!<"
  eventi.possibili <- colnames(dataLoaderOBJ$MMatrix)
  
  for( ID in names(dataLoaderOBJ$pat.process) ) {
    sequenza <- c("BEGIN",dataLoaderOBJ$pat.process[[ID]][,dataLoaderOBJ$csv.EVENTName],"END")
    # browser()
    for( ct in seq(1,length(sequenza)) ) {
      if(ct < (length(sequenza)) ){
        # if(sequenza[ ct ]=="Death") browser()
        EFPT.neverAfter[sequenza[ ct ] ,  sequenza[ (ct+1):length(sequenza)] ] <- "";
        mancanti <- eventi.possibili[ !(eventi.possibili %in% sequenza[ (ct+1):length(sequenza)]) ]
        EFPT.alwaysAfter[sequenza[ ct ] ,  mancanti ] <- "";
      }
      # if(ct >= 2){
      #   EFPT.neverBefore[sequenza[ ct ] ,  sequenza[ 1:(ct-1)] ] <- "";
      # }
    }
  }
  # EFPT.neverBefore[ "BEGIN" ,  ] <- "!<";
  EFPT.alwaysAfter[ "END" ,  ] <- "";
  return( list( "EFPT.hasNeverAfter" = EFPT.neverAfter,
                "EFPT.hasAlwaysAfter" = EFPT.alwaysAfter
                # "EFPT.neverBefore" = EFPT.neverBefore
                )
          )
}
  
#' calcola FT
#' 
#' @description  Funzione per calcolare la FM
#' @export
calcolaFootPrintTable.pat.process <- function( dataLoaderOBJ , skip.less.than = 0, threshold.perc = 0.00000001 ) {
  
  # rbind.data.frame <- do.call(rbind.data.frame, dataLoaderOBJ$pat.process )

  FPT<-dataLoaderOBJ$MMatrix;  FPT[,]<-"#"
  FPT.numbers.R<-dataLoaderOBJ$MMatrix;  FPT.numbers.R[,]<-0
  
  for( ID in names(dataLoaderOBJ$pat.process) ) {
    sequenza <- c("BEGIN",dataLoaderOBJ$pat.process[[ID]][,dataLoaderOBJ$csv.EVENTName],"END")
    # browser()
    for( ct in seq(1,length(sequenza)-1) ) {
      num <- FPT.numbers.R[ sequenza[ ct ] , sequenza[ (ct+1) ] ]
      FPT.numbers.R[ sequenza[ ct ] , sequenza[ (ct+1) ] ] <- num + 1
    }
  }
  
  FPT.numbers.R[  which(FPT.numbers.R < skip.less.than,arr.ind = T) ] <- 0
  
  ooo <- FPT.numbers.R
  for( riga in seq(1,nrow(ooo)) ) {
    for( colonna in seq(1,ncol(ooo)) ) {
      ooo[ riga, colonna ] <- abs(FPT.numbers.R[ riga, colonna ] )/abs(FPT.numbers.R[ riga, colonna ]+FPT.numbers.R[ colonna, riga ] )
    }
  }
  FPT.numbers.R.perc <- ooo
  FPT.numbers.R.perc[ which(is.nan(FPT.numbers.R.perc),arr.ind = T) ] <- 0
  
  
  for( riga in rownames(FPT.numbers.R.perc) ) {
    for( colonna in colnames(FPT.numbers.R.perc) ) {
      if( FPT.numbers.R.perc[ riga, colonna ] >= threshold.perc & 
          FPT.numbers.R.perc[ colonna, riga ] <= threshold.perc ) {
        FPT[ riga, colonna ] <- "->"
        FPT[ colonna, riga ] <- "<-"
      }
      if( FPT.numbers.R.perc[ riga, colonna ] >= threshold.perc & 
          FPT.numbers.R.perc[ colonna, riga ] <= threshold.perc ) {
        FPT[ riga, colonna ] <- "->"
        FPT[ colonna, riga ] <- "<-"
      }
      if( FPT.numbers.R.perc[ riga, colonna ] >= threshold.perc & 
          FPT.numbers.R.perc[ colonna, riga ] >= threshold.perc ) {
        FPT[ riga, colonna ] <- "||"
        FPT[ colonna, riga ] <- "||"
      }
      
    }
  }
  return( list("FPT"=FPT, "FPT.numbers.R"=FPT.numbers.R, "FPT.numbers.R.perc"=FPT.numbers.R.perc))
}

# -----------------------------------------------------------------------
# funzione plotPatientReplayedTimelineFunction
# -----------------------------------------------------------------------
#' Some useful tools new version
#' 
#' @description  A class which provide some tools. pMineR intarnal use only. wow
plotPatientReplayedTimelineFunction<-function( list.computation.matrix , patientID,
                                               text.cex=.7, y.intra.gap = 40, x.offset = 100,
                                               thickness=5 , 
                                               bar.border = "Navy",bar.volume = "lightsteelblue1",
                                               text.date.cex =.6) {
  
  date.notevoli <-c()
  durate.notevoli <- c()
  matrice <- list.computation.matrix$list.computation.matrix$stati.timeline[[patientID]]
  tempo.max <- max(  as.numeric(matrice[,4])  )
  numero.stati <- length(unique(matrice[,1]))
  arr.stati <- c()
  for( tmp in 1:length(matrice[,1])) {
    if( !(matrice[tmp,1] %in%arr.stati)) { arr.stati <- c(arr.stati,matrice[tmp,1]) }
  }
  # browser()
  par(mar=c(2,0,2,0)+0)
  plot( x=c(), y=c(), 
        xlim = c(0,tempo.max + x.offset+ 15) , 
        ylim=c(0,(numero.stati+1)*y.intra.gap ), 
        bty='n',axes = FALSE, xlab='', ylab='' )
  
  lista.boxes<- list()
  lista.points<- list()
  lista.date<-list()
  
  for( index in seq(1,length(arr.stati) )) {
    ypos.line <- (numero.stati+1)*y.intra.gap - index * y.intra.gap
    stato <- arr.stati[ index ]
    # text(x = 0,y = ypos.line,labels = stato, cex = text.cex, pos = 4)
    
    # lista.date[[length(lista.date)+1]] <- list("x"=c(x.offset,tempo.max+x.offset), "y"=c( ypos.line,ypos.line ))
    
    sub.matrice <- matrice[ which(matrice[ ,1]==stato )  ,]
    numero.righe.sub.matrice <- length(sub.matrice)/4
    # Se e' almeno una matrice (se ho almeno due rilevazioni)
    if(numero.righe.sub.matrice>1) {
      l.from <- NA
      l.to <- NA
      for( i in seq(1,numero.righe.sub.matrice )) {
        if(sub.matrice[i,2]=="begin") { 
          l.from <- as.numeric(sub.matrice[i,4]) 
          durate.notevoli <- c(durate.notevoli, l.from )
          date.notevoli <- c(date.notevoli, sub.matrice[i,3] )
          lista.date[[length(lista.date)+1]] <- list("x"=c(x.offset + l.from ,x.offset + l.from), "y"=c( -5, (numero.stati+1)*y.intra.gap +5),"label.data"=sub.matrice[i,3],"label.durata"=sub.matrice[i,4])          
        }
        if(sub.matrice[i,2]=="end") {
          l.to <- as.numeric(sub.matrice[i,4] )
          lista.date[[length(lista.date)+1]] <- list("x"=c(x.offset + l.to ,x.offset + l.to), "y"=c( -5, (numero.stati+1)*y.intra.gap +5),"label.data"=sub.matrice[i,3],"label.durata"=sub.matrice[i,4])          
          lista.boxes[[length(lista.boxes)+1]]<-list( "x"=c( l.from ,l.to, l.to, l.from, l.from ) + x.offset, "y"=c( -thickness, -thickness, thickness, thickness , -thickness)+ypos.line )
          durate.notevoli <- c(durate.notevoli, l.to )
          date.notevoli <- c(date.notevoli, sub.matrice[i,3]  )
        }
      }
    }
    # Se c'e' solo una riga!
    if(numero.righe.sub.matrice==1) {
      l.pos <- as.numeric(sub.matrice[4] )
      durate.notevoli <- c(durate.notevoli, l.pos )
      date.notevoli <- c(date.notevoli, sub.matrice[3]  )   
      lista.date[[length(lista.date)+1]] <- list("x"=c(x.offset + l.pos ,x.offset + l.pos), "y"=c( -5, (numero.stati+1)*y.intra.gap +5),"label.data"=sub.matrice[3], "label.durata"=sub.matrice[4])
      
      # Se e' un END 
      if(sub.matrice[2]=="end" |  as.numeric(sub.matrice[4])==tempo.max ) {
        lista.points[[length(lista.points)+1]]<-list("x"=l.pos + x.offset,"y"=ypos.line)
      }
      # Se e' un BEGIN
      if(sub.matrice[2]=="begin" & as.numeric(sub.matrice[4])!=tempo.max) {
        l.from <- l.pos
        l.to <- tempo.max
        lista.boxes[[length(lista.boxes)+1]]<-list( "x"=c( l.from ,l.to, l.to, l.from, l.from ) + x.offset, "y"=c( -thickness, -thickness, thickness, thickness , -thickness)+ypos.line )
      }
    }    
    
  }
  
  # plotta le verticali delle date
  number <- 1
  old.x <- c()
  for(i in seq(1, length(lista.date))) {
    if(! (lista.date[[i]]$x[1]  %in% old.x) ) {
      number <- number + 1 
      # points(x =lista.date[[i]]$x, y = lista.date[[i]]$y , type='l', col="grey", lty = 4 )
      points(x =lista.date[[i]]$x, y = lista.date[[i]]$y -15, type='l', col="grey", lty = 4 )
      text(x = lista.date[[i]]$x , y = lista.date[[i]]$y[1] + (number * 10)-5, labels = str_replace_all(string = lista.date[[i]]$label.data,pattern = " ",replacement = "\n"), cex = text.date.cex, col='black')
      text(x = lista.date[[i]]$x , y = (numero.stati+1)*y.intra.gap + (number * 10) -25, labels = as.integer(as.numeric(lista.date[[i]]$label.durata)), cex = text.date.cex, col='black')
      if(number >= 3) number <- 0
      old.x <- c(old.x, lista.date[[i]]$x[1] )
    }
  }
  # plotta gli assi degli stati
  for( index in seq(1,length(arr.stati) )) {
    points( x = c(x.offset,x.offset+tempo.max), 
            y = c( (numero.stati+1)*y.intra.gap - index * y.intra.gap, (numero.stati+1)*y.intra.gap - index * y.intra.gap),
            type='l' , col= "grey")     
  }
  # plotta i GANTT
  for(i in seq(1, length(lista.points))) {
    points( x = lista.points[[i]]$x, 
            y = lista.points[[i]]$y,
            pch=13 , col= bar.border)     
    # points(x =lista.date[[i]]$x, y = lista.date[[i]]$y , type='l', col="grey", lty = 4 )
  }  
  for(i in seq(1, length(lista.boxes))) {
    points( x = lista.boxes[[i]]$x, 
            y = lista.boxes[[i]]$y,
            type='l' , col= bar.border)
    polygon( x = lista.boxes[[i]]$x, 
             y = lista.boxes[[i]]$y,
             col= bar.volume) 
  }    
  
  for( index in seq(1,length(arr.stati) )) { 
    ypos.line <- (numero.stati+1)*y.intra.gap - index * y.intra.gap
    stato <- arr.stati[ index ]
    text(x = 0,y = ypos.line,labels = stato, cex = text.cex, pos = 4)
  }
  # list.computation.matrix
} 

#' A function to plot nice timeline
#'
#' @description  wow
#' @export
plotCohortTimeline <- function(  objDL.obj, lst.to.join=list() , lst.spike=list() , UM="days", y.grid = NA,
                                   ordered = FALSE, decreasing = TRUE , y.grid.col = "lightgrey", x.grid.col = "lightgrey") {
  
  
  # lst.to.join <- list( "Covid" = list( "from" ="Covid_BEGIN","to" = "Covid_END", "col"="red"), 
  #                      "Rehabilitation" = list( "from" = "Rehabilitation_BEGIN" , "to" = "Rehabilitation_END", col = "green"), 
  #                      "SubIntensive" = list( "from" = "SubIntensive_BEGIN" , "to" = "SubIntensive_END", col = "orange") 
  # )
  # lst.spike <- list( "Tested_Positive"=list( "col" = "blue", "lwd" = 2, "pch"= 8),
  #                    "Discharge"=list( "col" = "blue", "lwd" = 1, "pch"= 25)
  # )
  # plot.cohort.timeline(objDL.obj = objDL.new.export,lst.to.join = lst.to.join, lst.spike = lst.spike,y.grid = 90,
  #                      ordered = TRUE,decreasing = TRUE)
  
  
  
  objDL.new.export <- objDL.obj
  y.step.each <- y.grid
  
  if( UM == "mins") conversion <- 1
  if( UM == "hours") conversion <- 60
  if( UM == "days") conversion <- 24*60
  if( UM == "weeks") conversion <- 24*60*7
  
  
  minThickness <- 10
  internalMargin <- 2
  numberOfPatient <- length(objDL.new.export$pat.process)
  biggestTime <- max(do.call(rbind,objDL.new.export$pat.process)$pMineR.deltaDate)/conversion
  if(is.na(y.step.each)) y.step.each <- as.integer(biggestTime/10)
  
  max.time.per.patient <- unlist(lapply(names(objDL.new.export$pat.process), function(ID) { max(objDL.new.export$pat.process[[ID]]$pMineR.deltaDate) } ))
  
  if(ordered == TRUE) {
    ordered.patients <- order(max.time.per.patient,decreasing = !decreasing)
  } else { 
    ordered.patients <- 1:length(objDL.new.export$pat.process)
  }
  
  arr.posizioni <- seq( 0, (numberOfPatient+1)* minThickness, by= ((numberOfPatient+1)* minThickness) / numberOfPatient)
  arr.posizioni <- arr.posizioni[1:(length(arr.posizioni)-1)] + (minThickness/2)
  
  plot(0,0,xlim=c(0,biggestTime),ylim=c(0,(numberOfPatient+1)* minThickness ), axes=FALSE,xlab=UM, ylab="" , col="white")
  abline( v = seq(0,biggestTime,by = y.step.each) ,lty = 2, col=y.grid.col )
  
  tmp.arr.from <- unlist(lapply(names(lst.to.join),function(i) { lst.to.join[[i]]$from } ))
  tmp.arr.to <- unlist(lapply(names(lst.to.join),function(i) { lst.to.join[[i]]$to } ))
  tmp.arr.col <- unlist(lapply(names(lst.to.join),function(i) { lst.to.join[[i]]$col } ))
  arr.y.to.put.labels <- c()
  cursore <- 1
  tmp.1 <- lapply( ordered.patients, function( no.patient ) {
    y.start <- (cursore * minThickness)+internalMargin ; y.stop <- ((cursore+1)*minThickness)-internalMargin
    y.median <- (y.start + (y.stop-y.start)/2) 
    arr.y.to.put.labels <<- c( arr.y.to.put.labels , y.median )
    abline( h = y.median ,lty = 2, col=x.grid.col)
    arr.righe.2.skip <- c()
    tmp.2 <- lapply( 1:nrow(objDL.new.export$pat.process[[no.patient]]), function(riga) {

      if( !riga %in% arr.righe.2.skip) {
        arr.tutto <- objDL.new.export$pat.process[[no.patient]][[ objDL.new.export$csv.EVENTName ]]
        current.event <- objDL.new.export$pat.process[[no.patient]][[ objDL.new.export$csv.EVENTName ]][riga]
        browser()
        quale <- which( tmp.arr.from == current.event )
        if( length(quale)>0  ){
          da.trovare <- tmp.arr.to[quale]
          quale.to <- which(arr.tutto[(riga+1):length(arr.tutto)] == da.trovare)[1]
          quale.to <- quale.to + riga
          arr.righe.2.skip <- c( arr.righe.2.skip , quale.to )

          x.from <- (objDL.new.export$pat.process[[no.patient]]$pMineR.deltaDate[riga]/conversion)
          x.to <- (objDL.new.export$pat.process[[no.patient]]$pMineR.deltaDate[quale.to]/conversion)
          
          polygon( c(x.from , x.to, x.to , x.from , x.from  ),
                   c(y.start  , y.start ,y.stop , y.stop , y.start), lwd=2, 
                   col=tmp.arr.col[quale], border = NA)
          
        } else {
          if( current.event %in% names(lst.spike)) {
            
            points( (objDL.new.export$pat.process[[no.patient]]$pMineR.deltaDate[riga]/conversion), (y.start + (y.stop - y.start)/2),
                    col = lst.spike[[ current.event ]]$col, 
                    pch = lst.spike[[ current.event ]]$pch,
                    lwd = lst.spike[[ current.event ]]$lwd,
            )    
          }
        }      
      }
    } )
    cursore <<- cursore + 1
  })
  
  axis(2, arr.y.to.put.labels, labels=names(objDL.new.export$pat.process),las=2)
  axis(1, seq(0,biggestTime,by = y.step.each), labels=seq(0,biggestTime,by = y.step.each),las=1)
}

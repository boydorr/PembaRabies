#
## Function to distribute individuals to cells based on a probability map
#________________________________________


PopMap<- function(PembaUTM,PembaWard,init=0,probMap,wardPops,missedVill=NA,missedVillCells=NA,villageIDs){
  
  ## Set up column to hold initial population map
  PembaUTM$popMap <- init
  
  if(!is.null(missedVill)){
    if(is.na(missedVill)|is.na(missedVillCells)){
      
      ##Find villages that don't have an assigned cell
      missedVill <- which(!is.element(1:length(wardPops),unique(PembaUTM$WardID)))
      
      ##Find which grid cell each of these has the greatest degree of overlap with
      centroids <- SpatialPoints(coordinates(PembaWard[missedVill,]),proj4string=PembaUTM@proj4string)
      missedVillCells <- over(centroids,PembaUTM)$cellID
      
      ## If some missed villages haven't got a cell as they're on the coast, assign these to the closest cell
      if(length(which(is.na(missedVillCells)))>0){
        for (i in which(is.na(missedVillCells))){
          missedVillCells[i] <- which.min(gDistance(centroids[i,], PembaUTM, byid=TRUE))
        } 
      }
      
    } 
  }
  
    
  ## for each village
  for(i in 1:nrow(PembaWard@data)){
    
    if(is.element(i,missedVill)){
      PembaUTM$popMap[missedVillCells[which(missedVill==i)]] <- PembaUTM$popMap[missedVillCells[which(missedVill==i)]] + wardPops[i]
      
    }else {
      
      ## cells and individuals in village
      n_inds<-wardPops[i]
      if(sign(n_inds)==-1){
        cells<-which(villageIDs==i & PembaUTM$popMap>0)
        probs<-rep(1,length(cells))
      }else{
        cells<-which(villageIDs==i)
        probs <- probMap[which(!is.na(probMap[]))[cells]]
        if(sum(probs)==0){probs<-probs+1}
      }
      
      popChanges <- rep(0,length(cells))
      if(sign(n_inds)==-1){
        popChangesTable <- table(sample(rep(1:length(cells),times=PembaUTM$popMap[cells]),abs(n_inds)))
        popChanges[as.numeric(names(popChangesTable))] <- as.numeric(popChangesTable)
        
      }else if(sign(n_inds)==1){
        popChanges <- rmultinom(1,abs(n_inds),probs)
      }    
      
      
      ##distribute among available cells
      PembaUTM$popMap[cells] <- PembaUTM$popMap[cells] + sign(n_inds)*popChanges
      
    }
    
  }
  
  return(PembaUTM$popMap)
  
}



## Function for generating a randomized, incomplete block experimental design
### Tyler Tiede. University of Minnesota. May 22, 2014


# There are three different scenerios:
### 1) user-defined: the user enters the number of field rows and columns, as well as the number of rows and columns per block
### 2) semi-naive: the user only knows the number of field rows and the number of rows per blk
### 3) naive: the user only knows the number of field rows

#### There is not an option for starting with the number of columns because the case in
#### field research is typically that a researcher knows how long their plots are and how
#### wide the field is, w/ the length of the field available to plant varying

# This is a RETEST
# List of checks starting with the primary; assumes you have a primary check
design.dma<- function(enviro=format(Sys.Date(), "%x"), entries= NULL, nEntries= NULL, chk.names= NULL, nSecChk= NULL, nFieldRows= NULL, nFieldCols= NULL, nRowsPerBlk=NULL, nColsPerBlk=NULL, nChksPerBlk=2,  plot.start=1001, maxPerChks=0.12, fillWithChk=T, minRowBlkDim=2, minColBlkDim=3){
  
  ## Define user functions
  #JL Find all integers that are perfect divisors of x
  reduce <- function(x){
    keep <- c()
    for(i in 1:x){
      if(x %% i == 0) keep <- c(keep, i)   
    }
    return(keep)
  }
  #JL Same as reduce but exclude 1 and x.
  #JL If you are just testing for whether a number is prime, only go to floor(sqrt(x))
  #JL If function is to test whether prime or not, it should just return TRUE / FALSE
  prime <- function(x){
    keep <- c()
    for(i in 2:(x-1)){
      if(x %% i == 0) keep <- c(keep, i)   
    }
    return(keep)
  }
  
  
  ## QC of function inputs
  if(is.null(entries) & is.null(nEntries)) stop("Must provide an entry list (entries=) OR the number of entries desired (nEntries=).")
  
  if(is.null(chk.names) & is.null(nSecChk)) stop("Must provide a list of check names (chk.names=) with the primary check listed first\n OR the number of SECONDARY checks desired (nSecChk=).")
  
  if(is.null(nFieldRows)) stop("Provide the number of rows (sometimes called ranges or beds) (nFieldRows=).")
  
  if(is.null(prime(nFieldRows))) warning("nFieldRows is a prime number, thus rows will not be blocked.")
  ## Typically nFieldCols > nFieldRows so it is OK if you don't block rows as long as you have a sufficient number of columns to generate the number of blocks needed
  if(!is.null(nFieldCols)) if(is.null(prime(nFieldCols)) & (nFieldCols>nFieldRows)) stop("nFieldCols is a prime number and the columns cannot be blocked. Either change the input for/nnFieldCols to a non-prime number or increase the number of nFieldRows")
  
  ## Develop other non-input functions parameters
  if(is.null(entries)) entries <- as.matrix(paste("entry", 1:nEntries, sep="_"))
  
  if(is.null(nEntries)){
    entries <- as.matrix(entries) # If user input of entries was a list, will convert it to matrix
    nEntries <- nrow(entries)
  }
  
  # If the size of the field isn't large enough to accomdate all the entries + primary checks then it'll error
  #JL I think this should take the secondary checks into account also
  if(!is.null(nFieldCols)){ 
    if(((nFieldCols / nColsPerBlk) * (nFieldRows / nRowsPerBlk) + nEntries) > nFieldCols*nFieldRows){
      stop("The minimum number of plots has not been met by the given row/column dimensions.")
    } 
  }

  ## There is not a semi-naive procedure for a user input of nFieldRows, nRowsPerBlk, and nFieldCols... so this will
  ## provide an acceptable number of columns to be in a blk and use Scenerio 1
  if(!is.null(nFieldRows) & !is.null(nRowsPerBlk) & !is.null(nFieldCols) & is.null(nColsPerBlk)){
    nColsPerBlk <- minColBlkDim
    while(nFieldCols %% nColsPerBlk != 0) nColsPerBlk <- nColsPerBlk + 1
  }
  
  if(is.null(chk.names)){
    sec.chks <- as.character(2:(nSecChk+1)) ## Do need as seperate object for later on in function
    chk.names <- paste("chk", c(1,sec.chks), sep="") ## All generic check names
    rm(sec.chks) #JL I don't understand: if you need later, why do you remove here?
  }
  
  if(is.null(nSecChk)){
    sec.chks <- chk.names[-1] ## The primary check must be listed first in the function input
    nSecChk <- length(sec.chks)
    rm(sec.chks)
  }
  
  
  ##################################################################
  ######## Determining non-user submitted field paramaters #########
  ##################################################################
  
  ## Scenerio 1 (user-defined): User provides nFieldRows, nRowsPerBlk, nFldCols, nColsPerBlk
  
  if(!is.null(nFieldRows) & !is.null(nRowsPerBlk) & !is.null(nFieldCols) & !is.null(nColsPerBlk)){ ## If nFieldCols is defined
    
    if(nFieldRows %% nRowsPerBlk  != 0) stop("nFieldRows is not evenly divisible by nRowsPerBlk.")
    if(nFieldCols %% nColsPerBlk != 0) stop("nFieldCols is not evenly divisible by nColsPerBlk") # 5 cols/blk is the standard for the MADII design
    
    nBlkRows.tmp1 <- nFieldRows / nRowsPerBlk # 3 rows/block is the standard for the MADII design
    
    nBlkCols.tmp1 <- nFieldCols / nColsPerBlk # 5 cols/blk is the standard for the MADII design
    
    nBlks.tmp1 <- nBlkCols.tmp1 * nBlkRows.tmp1
    
    exp.size.tmp1 <- (nBlkRows.tmp1*nRowsPerBlk)*(nBlkCols.tmp1*nColsPerBlk)
    
    totEntries.tmp1 <- nEntries + (nBlks.tmp1*nChksPerBlk)
    
    ### Redefine tmp paramaters to final field paramaters
    nBlkRows <- nBlkRows.tmp1
    nBlkCols <- nBlkCols.tmp1
    nBlks <- nBlkRows.tmp1*nBlkCols.tmp1
     
  } #End of scenerio 1 loop
  
  ## Scenerio 2 (semi-naive): User provides nFieldRows and nRowsPerBlk
  if(!is.null(nFieldRows) & !is.null(nRowsPerBlk) & is.null(nFieldCols) & is.null(nColsPerBlk)){
    
    if(nFieldRows %% nRowsPerBlk  != 0) stop("nFieldRows is not evenly divisible by nRowsPerBlk.")
    
    nBlkRows.tmp2 <- nFieldRows / nRowsPerBlk # 3 rows/block is the standard for the MADII design
    
    
    finalPerChks <- 1
    startPerChks <- 0.10
    
    while(finalPerChks > maxPerChks){
      
      nChks.tmp2 <- ceiling((startPerChks * nEntries) / (1-startPerChks))
      
      totEntries.tmp2 <- nEntries + nChks.tmp2
      
      nFieldCols.tmp2 <- ceiling(totEntries.tmp2/nFieldRows)
      
      nBlks.tmp2 <- ceiling(nChks.tmp2/nChksPerBlk)
      
      ### Secondary checks need to be represented at least twice in the experiment
      while(nBlks.tmp2 < nSecChk*2) nBlks.tmp2 <- nBlks.tmp2
      
      nBlkCols.tmp2 <- ceiling(nBlks.tmp2/nBlkRows.tmp2)
      
      nColsPerBlk.tmp2 <- floor(nFieldCols.tmp2/nBlkCols.tmp2)
      
      exp.size.tmp2 <- (nBlkRows.tmp2*nRowsPerBlk)*(nBlkCols.tmp2*nColsPerBlk.tmp2)
      
      while(exp.size.tmp2 < totEntries.tmp2){
        nFieldCols.tmp2 <- nFieldCols.tmp2 + 1
        nColsPerBlk.tmp2 <- floor(nFieldCols.tmp2/nBlkCols.tmp2)
        exp.size.tmp2 <- (nBlkRows.tmp2*nRowsPerBlk)*(nBlkCols.tmp2*nColsPerBlk.tmp2) 
      }
      
      while(nFieldCols.tmp2 %% nColsPerBlk.tmp2 != 0 ) nFieldCols.tmp2 <- nFieldCols.tmp2 + 1
      nColsPerBlk.tmp2 <- floor(nFieldCols.tmp2/nBlkCols.tmp2)
      exp.size.tmp2 <- (nBlkRows.tmp2*nRowsPerBlk)*(nBlkCols.tmp2*nColsPerBlk.tmp2)
      
      if(fillWithChk) finalPerChks <- (exp.size.tmp2 - nEntries) / (exp.size.tmp2)
      if(!fillWithChk) finalPerChks <- nChks.tmp2/totEntries.tmp2
      startPerChks <- startPerChks - 0.0025
      
      
    }
    
    ### Redefine tmp paramaters to final field paramaters
    nFieldCols <- nFieldCols.tmp2
    nColsPerBlk <- nColsPerBlk.tmp2
    nBlkRows <- nBlkRows.tmp2
    nBlkCols <- nBlkCols.tmp2
    nBlks <- nBlkRows.tmp2*nBlkCols.tmp2
     
  } #End of scenerio 2 loop
  
  
  ## Scenerio 3 (naive): User provides ONLY nFieldRows
  if(!is.null(nFieldRows) & is.null(nRowsPerBlk) & is.null(nFieldCols) & is.null(nColsPerBlk)){
    
    ### Calculate starting (non-optimized paramaters)
    finalPerChks <- 1
    startPerChks <- 0.1
    
    while(finalPerChks >= maxPerChks){
      nChks.tmp3 <- ceiling((startPerChks * nEntries) / (1-startPerChks))
      
      totEntries.tmp3 <- nEntries + nChks.tmp3
      
      nFieldCols.tmp3 <- ceiling(totEntries.tmp3/nFieldRows)
      
      checkIfPrime <- prime(nFieldCols.tmp3)
      if(is.null(checkIfPrime)){
        nFieldCols.tmp3 <- nFieldCols.tmp3 + 1
        checkIfPrime <- prime(nFieldCols.tmp3)
      } ; rm(checkIfPrime)
      
      nBlks.tmp3 <- ceiling(nChks.tmp3/nChksPerBlk)
      
      #JL I think you are letting the tail wag the dog here
      #JL I would return an error saying that fewer secondary checks need to be given
      while(nBlks.tmp3 < nSecChk*2) nBlks.tmp3 <- nBlks.tmp3 + 1
      
      nBlks.min <- nBlks.tmp3
      
      nBlkRows.pool <- reduce(nFieldRows)[-c(1:(minRowBlkDim-1))]
      nBlkCols.pool <- reduce(nFieldCols.tmp3)[-c(1:(minColBlkDim-1))]
      
      rcList <- list(NULL) ; diffList <- c() ; count <-1
      for(r in nBlkRows.pool){
        for(c in nBlkCols.pool){
          nBlks.tmp <- r*c # nBlks.tmp different than nBlks.tmp3, just used within function
          rcList[[count]] <- c(r,c)
          diffList <- c(diffList, nBlks.tmp - nBlks.tmp3)
          count <- count + 1
        } 
      } ; rm(nBlks.tmp)
      
      pick <- which(diffList == min(diffList[diffList>=0]))
      
      diffList2 <- c()
      if(length(pick) == 1){
        nBlkRows.tmp3 <- rcList[[pick]][1]
        nBlkCols.tmp3 <- rcList[[pick]][2]
      } else{
        reduced.rcList <- rcList[pick]
        for(i in 1:length(reduced.rcList)){
          r <- reduced.rcList[[i]][1]
          c <- reduced.rcList[[i]][2]
          diffList2 <- c(diffList2, c-r)
        }
        pick2 <- which(diffList2 == min(diffList2[diffList2 >= 0]))
        
        nBlkRows.tmp3 <- reduced.rcList[[pick2]][1]
        nBlkCols.tmp3 <- reduced.rcList[[pick2]][2]
      }
      
      nRowsPerBlk.tmp3 <- nFieldRows/nBlkRows.tmp3
      nColsPerBlk.tmp3 <- nFieldCols.tmp3/nBlkCols.tmp3
      
      nBlks.tmp3 <- nBlkRows.tmp3*nBlkCols.tmp3
      nChks.tmp3 <- nBlks.tmp3*nChksPerBlk
      
      exp.size.tmp3 <- nFieldCols.tmp3*nFieldRows
      totEntries.tmp3 <- nEntries + nChks.tmp3
      
      ### If the original parameters are not able to encompass all of the entries then add a column and redo (below code)
      while(exp.size.tmp3 < totEntries.tmp3){
        #nblks.tmp3 <- nBlks.min
        nFieldCols.tmp3 <- nFieldCols.tmp3 + 1
        
        checkIfPrime <- prime(nFieldCols.tmp3)
        if(is.null(checkIfPrime)){
          nFieldCols.tmp3 <- nFieldCols.tmp3 + 1
          checkIfPrime <- prime(nFieldCols.tmp3)
        } ; rm(checkIfPrime)
        
        nBlkRows.pool <- reduce(nFieldRows)[-c(1:(minRowBlkDim-1))]
        nBlkCols.pool <- reduce(nFieldCols.tmp3)[-c(1:(minColBlkDim-1))]
        
        rcList <- list(NULL) ; diffList <- c() ; count <-1
        for(r in nBlkRows.pool){
          for(c in nBlkCols.pool){
            nBlks.tmp <- r*c # nBlks.tmp different than nBlks.tmp3, just used within function
            rcList[[count]] <- c(r,c)
            diffList <- c(diffList, nBlks.tmp - nBlks.tmp3)
            count <- count+1
          } 
        }
        
        pick <- which(diffList == min(diffList[diffList>=0]))
        
        diffList2 <- c()
        if(length(pick) == 1){
          nBlkRows.tmp3 <- rcList[[pick]][1]
          nBlkCols.tmp3 <- rcList[[pick]][2]
        } else{
          reduced.rcList <- rcList[pick]
          for(i in 1:length(reduced.rcList)){
            r <- reduced.rcList[[i]][1]
            c <- reduced.rcList[[i]][2]
            diffList2 <- c(diffList2, c-r)
          }
          pick2 <- which(diffList2 == min(diffList2[diffList2 >= 0]))
          
          nBlkRows.tmp3 <- reduced.rcList[[pick2]][1]
          nBlkCols.tmp3 <- reduced.rcList[[pick2]][2]
        }
        
        nRowsPerBlk.tmp3 <- nFieldRows/nBlkRows.tmp3
        nColsPerBlk.tmp3 <- nFieldCols.tmp3/nBlkCols.tmp3
        
        nBlks.tmp3 <- nBlkRows.tmp3*nBlkCols.tmp3
        nChks.tmp3 <- nBlks.tmp3*nChksPerBlk
        
        exp.size.tmp3 <- nFieldCols.tmp3*nFieldRows
        totEntries.tmp3 <- nEntries + nChks.tmp3
        
      }
      
      if(fillWithChk) finalPerChks <- (exp.size.tmp3 - nEntries) / (exp.size.tmp3)
      if(!fillWithChk) finalPerChks <- nChks.tmp3/totEntries.tmp3
      startPerChks <- startPerChks - 0.0025
      
    }
    
    nFieldCols <- nFieldCols.tmp3
    nRowsPerBlk <- nRowsPerBlk.tmp3
    nColsPerBlk <- nColsPerBlk.tmp3
    nBlkRows <- nBlkRows.tmp3
    nBlkCols <- nBlkCols.tmp3
    nBlks <- nBlkRows.tmp3*nBlkCols.tmp3
    
  } #End of scenerio 3 loop
  
  
  ##################################################################
  ################### Finalize Paramaters ##########################
  ##################################################################
  
  ### Now that the field dimension are set
  exp.size <- nFieldRows*nFieldCols
  
  nPlotsPerBlk <- nRowsPerBlk * nColsPerBlk
  nChks <- nBlks*nChksPerBlk
    
  nSecChkPlots <- nBlks*(nChksPerBlk-1) 
  nSecChkReps <- floor(nSecChkPlots/nSecChk) # Number of times each secondary check will be observed in the field
  nSecChkPerBlk <- nChksPerBlk - 1 # Always have a primary check in block... the rest of checks are primary
    
  chk2.names <- chk.names[-1]
  chk2.nums <- 2:(length(chk2.names)+1)
  
  secChkPool <- rep(chk2.nums, times=nSecChkReps)
  secChkPool <- c(secChkPool, sample(chk2.nums, size=nSecChkPlots-length(secChkPool), replace=F)) 
  secChkPool <- secChkPool[order(randu[[1]][1:length(secChkPool)])]
   
  nFill <- exp.size - (nEntries + nBlks*nChksPerBlk) # Fill lines are empty plots at the end of the experiment
  
  if(!fillWithChk & nFill > 0){
    if(nFill >= nSecChk){ ## This part depends on user inputs so will not automatically adjust anything, just provide a warning
      warning("\nAn excessive number of fill plots are currently being used. It is recommended to:\n1) increase the nChksPerBlk, 2) decrease block size, or 3) decrease experiment size. \nYour design has still been created\n ")
    }
  }
    
  
  #################################################################
  ############### Build Field File ################################
  #################################################################
  
  fld.dgn <- as.data.frame(matrix(data=c((plot.start:(plot.start-1+exp.size)), rep(1:nFieldRows, each=nFieldCols), rep(c(1:nFieldCols, nFieldCols:1), length.out=exp.size), rep(NA, times=exp.size), rep(1:nBlkRows, each=(exp.size/nBlkRows)), rep(c(1:nBlkCols, nBlkCols:1), each=5, length.out=exp.size), rep(NA, times=2*exp.size)), nrow=exp.size, byrow=F))
  colnames(fld.dgn) <- c("Plot", "Row", "Col", "Blk", "Row.Blk", "Col.Blk", "Line.Code", "Entry")
  
  if(!fillWithChk & nFill>0){
    fld.dgn[(exp.size-nFill+1):exp.size, 7] <- "F"
  }
  
  blk.list <- 1:nBlks
  
  for(b in 1:nBlkRows){
    if((b %% 2 == 0)){
      blk.list[(1+nBlkCols*(b-1)):((nBlkCols*(b-1))+nBlkCols)] <- rev((1+nBlkCols*(b-1)):((nBlkCols*(b-1))+nBlkCols)) 
    } else{
      blk.list[(1+nBlkCols*(b-1)):((nBlkCols*(b-1))+nBlkCols)] <- ((1+nBlkCols*(b-1)):((nBlkCols*(b-1))+nBlkCols))
    }
  }
  
  ## Assign plots to blocks
  count <- 1
  for(b in 1:nBlkRows){
    for(c in 1:nBlkCols){
      blk <- blk.list[count]
      r1 <- (1+seq(0,1000,by=nRowsPerBlk))[b]
      r2 <- (seq(0,1000,by=nRowsPerBlk))[b+1]
      
      c1 <- (1+seq(0,250000,by=nColsPerBlk))[c]
      c2 <- (seq(0,250000,by=nColsPerBlk))[c+1]
      
      fld.dgn[which(fld.dgn$Row %in% r1:r2 & fld.dgn$Col %in% c1:c2), 4] <- blk
      
      count <- count+1
    } 
  }
  
  ## Assign primary checks first AND secondary checks
  
  for(b in 1:nBlks){
    blk <- fld.dgn[which(fld.dgn$Blk==b), ]
    blk <- blk[which(is.na(blk$Line.Code)),]
    row <- which(fld.dgn$Row == floor(mean(blk$Row)))
    col <- which(fld.dgn$Col == floor(mean(blk$Col)))
    fld.dgn[row[which(row  %in% col)], 7] <- 1 
  }
  
  for(bb in 1:nBlks){
    blk <- fld.dgn[which(fld.dgn$Blk==bb), ]
    blk2 <- blk[which(is.na(blk$Line.Code)),]
    samp <- sample(1:nrow(blk2), 1)
    row2 <- blk2[samp, "Row"]
    col2 <- blk2[samp, "Col"]
    fld.dgn[which(fld.dgn$Row==row2 & fld.dgn$Col==col2), 7] <- secChkPool[bb]   
  }
  
  #### Instead of fill, add secondary checks randomly to the experiment
  fld.dgn.save <- fld.dgn
  if(fillWithChk & (nFill > 0)){
    
    satisfied <- F
    
    while(!satisfied){
      fld.dgn <- fld.dgn.save
      tmp.min <- 10
      fld.dgn <- fld.dgn.save
      #fld.dgn[which(fld.dgn$Line.Code=="F"), "Line.Code"] <- NA
      chkFills <- sample(rep(chk2.nums, times=ceiling(nFill/length(chk2.nums))), nFill, replace=F)
      
      
      
      fld.plots.sel <- sample(as.numeric(rownames(fld.dgn[which(is.na(fld.dgn$Line.Code)),])), nFill)
      fld.dgn[fld.plots.sel, 7] <- chkFills  
      
      for(c in chk2.nums){
        tmp <- length(which(fld.dgn$Line.Code==c))
        if(tmp < tmp.min) tmp.min <- tmp
      }
      if(tmp.min >= 2) satisfied <- T
    }
  }

  ## Entries are coded 0
  fld.dgn[which(is.na(fld.dgn$Line.Code)), 7] <- 0
  
  options(warn=-1) ## If no fill plots then the lines below will throw an error
  ## Assign entry names and check names
  fld.dgn[which(fld.dgn$Line.Code == 0), 8] <- entries[order(sample(1:length(entries), length(entries)))]
  
  if(!fillWithChk & nFill>0){
    fld.dgn[which(fld.dgn$Line.Code == "F"), 8] <- "Fill"
  }
  fld.dgn[which(fld.dgn$Line.Code == "F"), 7] <- 0 # Change the fill lines back to experimental entries for purposes of MADIIadj
  options(warn=0) # Turn warnings back on
  
  for(c in 1:length(chk.names)){
    chk <- chk.names[c]
    fld.dgn[which(fld.dgn$Line.Code==c), 8] <- chk
  }
  
  ## Ensure that each secondary check is represented at least twice
  ## This will only be an issue under Scenerio 1
  tmp.min <- 100
  for(c in chk2.nums){
    tmp <- length(which(fld.dgn$Line.Code==c))
    if(tmp < tmp.min) tmp.min <- tmp
  }
  if(tmp.min < 2) warning("One or more of the secondary checks is represented only once in the experiment./nReduce the number of secondary checks.")
    
  rlzPerChks <- 1-(nrow(fld.dgn[which(fld.dgn$Line.Code=="0"),])/nrow(fld.dgn))
  
  ## Instead of re-working code, just re-order fld.dgn df before reading out
  fld.dgn.mod <- cbind(matrix(enviro, nrow=nFieldRows, ncol=1), fld.dgn) ; colnames(fld.dgn.mod)[c(1,8)] <- c("Enviro", "Check")
  fld.dgn.mod <- fld.dgn.mod[,c(1,2,9,8, 3:7)]
  fld.dgn$Line.Code
  
  design.ma<-list()
  
  design.ma[[1]]<-fld.dgn
  design.ma[[2]]<-fld.dgn.mod
  design.ma[[3]]<-nSecChk
  design.ma[[4]]<-chk.names
  design.ma[[5]]<-nFill
  design.ma[[6]]<-nBlkRows
  design.ma[[7]]<-nBlkCols
  design.ma[[8]]<-rlzPerChks
  
  return(design.ma)
  
} ## End of design.dma function






# The wikipedia entry on AR1 (Autoregressive_model) gives 
# var = \sigma^2_e / (1 - \phi^2) and autocovariance as cov = \phi^dist * \sigma^2_e / (1 - \phi^2)
# So \phi is the rate of decay of the error covariance with physical distance:
# 0 < \phi < 1.  \phi close to zero means little covariance from one plot to the next. 
# \phi close to one means error deviations change little from one plot to the next.
# This function parameters: field size in terms of numbers of rows and columns; plot size; \phi; \sigma^2_e
# You can supply a MADII design supplied by MADIIdgn instead of nRows and nCols
# plotSize is a two-vector with plotSize[1] the width (along the rows) and plotSize[2] the length (along the columns)
# ar1CovMat returns a covariance matrix that follows the autoregressive model with a variance of 1.
# That is, sigma^2_e = 1 - \phi^2 and so \sigma^2_e does not need to be supplied.
# if a MADII design was supplied (object created by MADIIdgn), matrix rows and columns 
# correspond to the plot order given in the design.
# if nRows and nCols were supplied, the plot order could be formed into a matrix by matrix(plotOrder, nRows, nCols)

ar1CovMat <- function(madII=NULL, nRows=NULL, nCols=NULL, plotSize=c(1, 2), phi=0.5){
  if (!is.null(madII)){
    rowVec <- madII$Row
    colVec <- madII$Col
    nRows <- max(rowVec)
    nCols <- max(colVec)
  } else{
    rowVec <- rep(1:nRows, times=nCols)
    colVec <- rep(1:nCols, each=nRows)    
  }
  # Calculate the distance matrix
  rowDistMat <- sapply(rowVec, function(rowNum) abs(rowNum - rowVec))
  colDistMat <- sapply(colVec, function(colNum) abs(colNum - colVec))
  distMat <- sqrt(plotSize[1]^2 * rowDistMat^2 + plotSize[2]^2 * colDistMat^2)
  return(phi^distMat)
}

# Function to determine the best decay of autocorrelation and predict entry effects
# The data.frame fieldObs should have (at least) two columns, Entry which should be the same as the madII Entry
# and phenoVal, which are the observations from the field
# The function tries \phi values from 0.1 to 0.9 (\phi is the rate of decay of the error covariance with distance)
ar1analysis <- function(madII, fieldObs, plotSize=c(1, 2)){
  maxLLik <- -Inf
  phiVals <- seq(from=0.1, to=0.9, by=0.1)
  for (phi in phiVals){
    regData <- list(phenoVal=fieldObs$phenoVal, Entry=fieldObs$Entry, R=ar1CovMat(madII, plotSize=plotSize, phi=phi))
    regOut <- regress(phenoVal ~ 1, ~ Entry + R, pos=c(TRUE, TRUE), identity=FALSE, data=regData)
    if (regOut$llik > maxLLik){
      bestReg <- list(regOut=regOut, phi=phi)
      maxLLik <- regOut$llik
    }
  }
  return(bestReg)
}

# This function does the same thing as ar1analysis but uses the R function optim.  I'm not sure what goes faster.
ar1optim <- function(madII, fieldObs, plotSize=c(1, 2)){
  toMinimize <- function(phi){
    regData <- list(phenoVal=fieldObs$phenoVal, Entry=fieldObs$Entry, R=ar1CovMat(madII, plotSize=plotSize, phi=phi))
    return(-regress(phenoVal ~ 1, ~ Entry + R, pos=c(TRUE, TRUE), identity=FALSE, data=regData)$llik)
  }
  optAR1 <- optim(par=0.5, fn=toMinimize, method="Brent", lower=0, upper=1)
  regData <- list(phenoVal=fieldObs$phenoVal, Entry=fieldObs$Entry, R=ar1CovMat(madII, plotSize=plotSize, phi=optAR1$par))
  return(list(regOut=regress(phenoVal ~ 1, ~ Entry + R, pos=c(TRUE, TRUE), identity=FALSE, data=regData), phi=optAR1$par))
}

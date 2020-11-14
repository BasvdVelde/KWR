#Similarity functions
sim.func <- function(x, y, algorithm = NULL) {
  Check <- c("NDP", "PC", "SC", "JI")
  if (is.null(algorithm)) {
    message("No 'algorithm' defined. Using 'NDP' (Normalized dot product) by default")
    algorithm <- "NDP"
  } 
  if (!algorithm %in% Check) message("Algorithm must be 'NDP' (Normalized dot product), 'PC' (Pearson's correlation), 'SC' (Spearman's correlation), or 'JI' (Jaccard index)")
  
  if(algorithm == "NDP") {
    as.vector((x %*% y)/(sqrt(sum(x^2)) * 
                         sqrt(sum(y^2))))
  } else if(algorithm == "PC") {
      as.vector((sum((x-mean(x))*(y-mean(y))))/
                (sqrt(sum((x-mean(x))^2)) *
                   sqrt(sum((y-mean(y))^2))))
  } else if(algorithm == "SC") {
      cor(x,y, method = "spearman")
      
  } else if(algorithm == "JI") {
      as.vector(length(intersect(x, y))/length(union(x, y)))
    }
} 

#Function that generates random data for the RawData frame to generate unrelated pairs
Randomize <- function(RawData) {
  ColTP <- RawData[,1:11]
  ColPC <- RawData[,12:22]
  
  shuffleTP <- ColTP[sample(nrow(ColTP)),]
  shufflePC <- ColPC[sample(nrow(ColPC)),]
  
  RawData <- cbind(shuffleTP, shufflePC)
}

#root mean squared error calculation
MSE <- function(x, y) {
  (sum(((x-y)^2)/length(x)))
}

#R2 
rsq <- function(x,y) {
  cor(x, y) ^ 2
}

#R2 with linear line (x +0)
R2L <- function(x,y) {
  yhat <- (1*x+0)
  RSS <- sum((y - yhat)^2)
  TSS <- sum((y - mean(y))^2)
  R2L <- 1-(RSS/TSS)
  R2L
}

#binning
binning <- function(spec.top, spec.bottom, xlim, b, t) { 
top_tmp <- data.frame(mz = spec.top[, 1], intensity = spec.top[,2]) #Creates a dataframe when the mz values are in the first column and the intensity in the second column
top_tmp$normalized <- ((top_tmp$intensity/max(top_tmp$intensity)) * 100) #creates a new column where the intensity is normalized to 100
top_tmp <- subset(top_tmp, top_tmp$mz >= xlim[1] & top_tmp$mz <= xlim[2]) #checks if the xlim set in the function is lower and higher than the selected mz values to ensure all values lie in the range
top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized) #creates a data.frame with only mz and normalized intensity - non-normalized intensity is discarded
top <- subset(top_plot, top_plot$intensity >= b) #keeps the normalised intensities higher than the set b value, which in this case was 10

bottom_tmp <- data.frame(mz = spec.bottom[, 1], intensity = spec.bottom[,2]) #Creates a dataframe when the mz values are in the first column and the intensity in the second column
bottom_tmp$normalized <- ((bottom_tmp$intensity/max(bottom_tmp$intensity)) * 100)  #creates a new column where the intensity is normalized to 100
bottom_tmp <- subset(bottom_tmp, bottom_tmp$mz >= xlim[1] & bottom_tmp$mz <= xlim[2]) #checks if the xlim set in the function is lower and higher than the selected mz values to ensure all values lie in the range
bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized) #creates a data.frame with only mz and normalized intensity - non-normalized intensity is discarded
bottom <- subset(bottom_plot, bottom_plot$intensity >= b) #keeps the normalised intensities higher than the set b value, which in this case was 10


for (i in 1:nrow(bottom)) top[, 1][bottom[, 1][i] >= top[, 1] - t & bottom[, 1][i] <= top[, 1] + t] <- bottom[,1][i]
alignment <- merge(top, bottom, by = 1, all = TRUE)
if (length(unique(alignment[, 1])) != length(alignment[,1])) 
  warning("the m/z tolerance is set too high")
alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0
names(alignment) <- c("mz", "intensity.top", 
                      "intensity.bottom")
alignment
} 

#collecting JSON file from MBID
MBID <- function(MBID, precursor = FALSE) {
  if (is.na(MBID))  
    return(data.frame())
  
  json_file <- paste("E:/KWR data/MBID", "/", MBID, ".JSON", sep = "")
  result <- fromJSON(file = json_file)
  peaks <- result$peaks
  MW <- result$precursor_mz
  MSintensity <- lapply(peaks, `[[`, 2)
  MSmz <- lapply(peaks, `[[`, 1)
  MSintensity <- as.data.frame(unlist(MSintensity))
  MSmz <- as.data.frame(unlist(MSmz))
  spec <- cbind(MSmz, MSintensity)
  colnames(spec) <- c("mz", "intensity")
  
  if (precursor == TRUE) { 
    
  return(MW) 
  
  } else 
    
  if (precursor == FALSE) { 
  
  return(spec)
  }
} 

allSame <- function(l)
{
  if (length(l) > 1)
  {
    if (all(is.na(l)))
      return(TRUE)
    if (any(is.na(l)))
      return(FALSE)
    
    return(all(l[[1]] == l))
  }
  
  return(TRUE)
}

classification <- function(x, y, classification = "Maccs") {
  if (classification == "Maccs") {
    RawData <- rbind(x, y)
    class1 <- RawData[RawData$Maccs >= 0.5,]
    class1$class <- 1
    class0 <- RawData[RawData$Maccs < 0.5,]
    class0$class <- 0
    RawData <- rbind(class1, class0)
    
  } else if (classification == "Known") {
    
    RawData1$class <- 1
    RawData2$class <- 0
    RawData <- rbind(RawData1, RawData2)
    
  }
}

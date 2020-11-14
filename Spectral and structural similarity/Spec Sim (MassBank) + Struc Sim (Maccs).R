################################################################################
#                             Collect data from MassBank                       #
################################################################################

source("utils.R")
RawData1 <- read.csv("Related_pairs.csv")
RawData2 <- read.csv("Random_pairs.csv")

RawData <- classification(RawData1, RawData2, classification = "Known")

library("rjson")


{ 
  pb   <- txtProgressBar(1, nrow(RawData), style=3)
  
  allResults <- mapply(RawData$ME, RawData$ME.1, RawData$SMILES.TP, RawData$SMILES.PC, RawData$class, seq_len(nrow(RawData)), FUN = function(MBIDTP, MBIDPC, SMILESTP, SMILESPC, class, i) {
    if (is.na(MBIDTP) || is.na(MBIDPC))
      return(NULL)
    
    setTxtProgressBar(pb, i)
    #Parent Compound
    #json_filePC <- sprintf("https://metabolomics-usi.ucsd.edu/json/?usi=mzspec:MASSBANK:%s",MBIDPC) #When extraction from online (needed when the JSON files are not op PC)
    json_filePC <- paste("MBID", "/", MBIDPC, ".JSON", sep = "")
    resultPC <- fromJSON(file = json_filePC)
    
    peaksPC <- resultPC$peaks
    MWPC <- resultPC$precursor_mz
    
    MSintensityPC <- lapply(peaksPC, `[[`, 2)
    MSmzPC <- lapply(peaksPC, `[[`, 1)
    
    MSintensityPC <- as.data.frame(unlist(MSintensityPC))
    MSmzPC <- as.data.frame(unlist(MSmzPC))
    
    spec.top <- cbind(MSmzPC, MSintensityPC)
    colnames(spec.top) <- c("mz", "intensity")
    
    #if (nrow(spec.top) == 0) browser()
    #if (nrow(spec.bottom) == 0) browser()
    
    #Transformation product
    #json_fileTP <- sprintf("https://metabolomics-usi.ucsd.edu/json/?usi=mzspec:MASSBANK:%s",MBIDTP) #When extraction from online (needed when the JSON files are not op PC)
    json_fileTP <- paste("MBID", "/", MBIDTP, ".JSON", sep = "")
    resultTP <- fromJSON(file = json_fileTP)
    
    peaksTP <- resultTP$peaks
    MWTP <- resultTP$precursor_mz
    
    MSintensityTP <- lapply(peaksTP, `[[`, 2)
    MSmzTP <- lapply(peaksTP, `[[`, 1)
    
    MSintensityTP <- as.data.frame(unlist(MSintensityTP))
    MSmzTP <- as.data.frame(unlist(MSmzTP))
    
    spec.bottom <- cbind(MSmzTP, MSintensityTP)
    colnames(spec.bottom) <- c("mz", "intensity")
    
    ################################################################################
    #                                Set parameters                                #
    ################################################################################
    
    b <-  0.5 #In the alignment it selects the minimum normalized intensity. When set to 0 all peaks are considered 
    
    t <- 0.005  #The mass error for alignment 
  
    c1 <- 0.0
    d1 <- 1.0
    
    c2 <- 0.0
    d2 <- 1.0
    
    c3 <- 0.0
    d3 <- 1.0
    
    xlim = c(50, max(MWPC, MWTP)+10)
    
    #Some default settings are:
    #NIST:        c=3 d=0.6
    #MASSBANK     c=2 d=0.5
    
    ################################################################################
    #                     Spectral similarity   (Approach 1)                       #
    ################################################################################
    
    #bin the two MS/MS spectra 
    alignment <- binning(spec.top, spec.bottom, xlim, b, t)
    
    alignment <- alignment[!findInterval(alignment$mz, c(MWPC -t, MWPC +t)) == 1,] #removes precursor PC
    alignment <- alignment[!findInterval(alignment$mz, c(MWTP -t, MWTP +t)) == 1,] #removes precursor TP
    
    mz <- alignment[, 1]
    x1 <- alignment[, 2]
    x2 <- alignment[, 3]
    
    u <- as.vector(mz^c1 * x1^d1)
    v <- as.vector(mz^c1 * x2^d1)
    
    #Normalized dot.product
    similarity_scoreDP <- sim.func(u, v, algorithm = "NDP")
    
    
    # if (allSame(u)) browser()
    #Pearson's correlation
    Similarity_scoreC <- if (length(u) < 3 || length(v) < 3) {
      sim.func(u, v, algorithm = "SC")
    }  else if (!allSame(u) && !allSame(v) && shapiro.test(u)$p.value >= 0.05 && shapiro.test(v)$p.value >= 0.05) {
      sim.func(u, v, algorithm = "PC") 
    } else { 
      sim.func(u, v, algorithm = "SC")
    }
    
    Similarity_scoreC <- abs(Similarity_scoreC)
    
    #Jaccard index
    jx1 <- alignment[(alignment[,2]>0),]$mz
    jx2 <- alignment[(alignment[,3]>0),]$mz
    
    Similarity_scoreJI <- sim.func(jx1, jx2, algorithm = "JI")
    
    #data for barplot
    data <- data.frame(algorithm=c("Normalized dot product", "Pearon's correlation", "Jaccard index"),  
                       approach1=c(similarity_scoreDP, Similarity_scoreC, Similarity_scoreJI))
    
    ################################################################################
    #              Spectral similarity fragment shift  (Approach 2)                #
    ################################################################################
    spec.topff <- spec.top[!findInterval(spec.top$mz, c(MWPC -t, MWPC +t)) == 1,] #removes precursor PC
    spec.bottomff <- spec.bottom[!findInterval(spec.bottom$mz, c(MWTP -t, MWTP +t)) == 1,] #removes precursor TP
    
    diff=(MWPC - MWTP)
    
    fragments.top <- spec.topff$mz 
    fragments.bottom <- (spec.bottomff$mz)+diff #orginal intensities but with mass difference to precusor peak
    
    spec.topf <- data.frame("mz" = fragments.top, "intensity" = spec.topff$intensity)
    spec.bottomf <- data.frame("mz" = fragments.bottom, "intensity" = spec.bottomff$intensity)
    
    #bin the two MS/MS spectra 
    alignmentf <- binning(spec.topf, spec.bottomf, xlim, b, t)
    
    mzf <- alignmentf[, 1]
    x1f <- alignmentf[, 2]
    x2f <- alignmentf[, 3]
    
    uf <- as.vector(mzf^c2 * x1f^d2)
    vf <- as.vector(mzf^c2 * x2f^d2)
    
    #Normalized dot.product
    similarity_scoreDPf <- sim.func(uf, vf, algorithm = "NDP")
    
    #Pearson's correlation
    Similarity_scoreCf <- if (length(u) < 3 || length(v) < 3) {
      sim.func(uf, vf, algorithm = "SC")
    }  else if (!allSame(u) && !allSame(v) && shapiro.test(u)$p.value >= 0.05 && shapiro.test(v)$p.value >= 0.05) {
      sim.func(uf, vf, algorithm = "PC") 
    } else { 
      sim.func(uf, vf, algorithm = "SC")
    }
    
    Similarity_scoreCf <- abs(Similarity_scoreCf)
    
    #Jaccard index
    jx1f <- alignmentf[(alignmentf[,2]>0),]$mz
    jx2f <- alignmentf[(alignmentf[,3]>0),]$mz
    
    Similarity_scoreJIf <- sim.func(jx1f, jx2f, algorithm = "JI")
    
    #data for barplot
    dataf <- data.frame(algorithm=c("Normalized dot product", "Pearon's correlation", "Jaccard index"),  
                        approach2=c(similarity_scoreDPf, Similarity_scoreCf, Similarity_scoreJIf))
    
    ################################################################################
    #             Spectral similarity fragment shift 2 (Approach 3)                #
    ################################################################################
    intersect.v <- as.vector(intersect(jx1, jx2)) #looks which peaks align
    
    alignmentff2 <- alignment
    
    #alignmentff2 <- alignmentff2[!findInterval(alignmentff2$mz, c(MWPC -t, MWPC +t)) == 1,] #removes precursor PC
    #alignmentff2 <- alignmentff2[!findInterval(alignmentff2$mz, c(MWTP -t, MWTP +t)) == 1,] #removes precursor TP
    
    adj.alignment <- alignmentff2[!alignmentff2$mz %in% c(intersect.v),] #alignment that contains the removed precursor ions and the aligned spectra are also removed
    
    mz.bottom <- adj.alignment$mz + diff #orginal intensities but with mass difference to precursor peak
    
    topf2 <- adj.alignment[,1:2]
    topf2 <- topf2[!(apply(topf2, 1, function(y) any(y == 0))),]
    bottomf2 <- cbind(data.frame(mz.bottom), data.frame(adj.alignment$intensity.bottom))
    bottomf2 <- bottomf2[!(apply(bottomf2, 1, function(y) any(y == 0))),]
    
    for (i in 1:nrow(bottomf2)) topf2[, 1][bottomf2[, 1][i] >= topf2[, 1] - t & bottomf2[, 1][i] <= topf2[, 1] + t] <- bottomf2[,1][i]
    alignmentf2 <- merge(topf2, bottomf2, by = 1, all = TRUE)
    if (length(unique(alignmentf2[, 1])) != length(alignmentf2[,1])) 
      warning("the m/z tolerance is set too high")
    alignmentf2[, c(2, 3)][is.na(alignmentf2[, c(2, 3)])] <- 0
    names(alignmentf2) <- c("mz", "intensity.top", 
                            "intensity.bottom")
    
    adj.alignment1 <- alignment[alignment$mz %in% c(intersect.v),] #this contains all the aligned spectra
    
    alignmentf2.unordered <- rbind(alignmentf2, adj.alignment1) #alignmentf2 has all the not aligned spectra with mass difference
    
    alignmentf3 <- alignmentf2.unordered[order(alignmentf2.unordered$mz),]
    
    mzf2 <- alignmentf3[, 1]
    x1f2 <- alignmentf3[, 2]
    x2f2 <- alignmentf3[, 3]
    
    uf2 <- as.vector(mzf2^c3 * x1f2^d3)
    vf2 <- as.vector(mzf2^c3 * x2f2^d3)
    
    #Normalized dot.product
    similarity_scoreDPf2 <- sim.func(uf2, vf2, algorithm = "NDP")
    
    #Pearson's correlation
    Similarity_scoreCf2 <- if (length(u) < 3 || length(v) < 3) {
      sim.func(uf2, vf2, algorithm = "SC")
    }  else if (!allSame(u) && !allSame(v) && shapiro.test(u)$p.value >= 0.05 && shapiro.test(v)$p.value >= 0.05) {
      sim.func(uf2, vf2, algorithm = "PC") 
    } else { 
      sim.func(uf2, vf2, algorithm = "SC")
    }
    
    Similarity_scoreCf2 <- abs(Similarity_scoreCf2)
    
    
    p1 <- if (length(u) < 3 || allSame(u)) { 
      p1 = 0.00 } 
    else { 
      p1 <- shapiro.test(u)$p.value
    }
    
    p2 <- if (length(v) < 3 || allSame(v)) { 
      p2 = 0.00 } 
    else { 
      p2 <- shapiro.test(v)$p.value
    }
    
    #Jaccard index
    jx1f2 <- alignmentf3[(alignmentf3[,2]>0),]$mz
    jx2f2 <- alignmentf3[(alignmentf3[,3]>0),]$mz
    
    Similarity_scoreJIf2 <- sim.func(jx1f2, jx2f2, algorithm = "JI")
    
    #data for barplot
    dataf2 <- data.frame(algorithm=c("Normalized dot product", "Pearon's correlation", "Jaccard index"),  
                         approach3=c(similarity_scoreDPf2, Similarity_scoreCf2, Similarity_scoreJIf2))
    
    ###############################################################################
    #                                    Tanimoto sim                              #
    ################################################################################
    SMILES <- cbind(SMILESTP, SMILESPC)
    library("rcdk")
    
    mols <- parse.smiles(SMILES, kekulise = TRUE)
    
    for (m in mols)
    {
      rcdk::do.typing(m)
      rcdk::do.aromaticity(m)
    }
    
    fps <- lapply(mols, get.fingerprint, type="maccs")
    
    fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')
    
    ################################################################################
    #                             dataframe for all results                        #
    ################################################################################
    Results <- cbind(data, "approach2" = dataf[,2], "approach3" = dataf2[,2], "Tanimoto_sim" = fp.sim[2,1], "class" = class, "p1" = p1, "p2" = p2)
    
    return(Results)
  }, SIMPLIFY = FALSE)
  
  
  
  names(allResults) <- paste0("name" = RawData$Transformation_Product, "&", RawData$Parent_compound)
  allResults <- allResults[!sapply(allResults, is.null)]
  
  allResultsDF <- do.call(rbind, allResults)
  
  ################################################################################
  #                                     Box plot                                 #
  ################################################################################
  allResultsNDP <-  allResultsDF[(allResultsDF$algorithm == "Normalized dot product"),]
  allResultsNDP$approach1[is.nan(allResultsNDP$approach1)] <- 0
  allResultsNDP$approach2[is.nan(allResultsNDP$approach2)] <- 0
  allResultsNDP$approach3[is.nan(allResultsNDP$approach3)] <- 0
  
  allResultsC <-  allResultsDF[(allResultsDF$algorithm == "Pearon's correlation"),]
  allResultsC$approach1[is.nan(allResultsC$approach1)] <- 0
  allResultsC$approach2[is.nan(allResultsC$approach2)] <- 0
  allResultsC$approach3[is.nan(allResultsC$approach3)] <- 0
  
  allResultsJI <-  allResultsDF[(allResultsDF$algorithm == "Jaccard index"),]
  allResultsJI$approach1[is.nan(allResultsJI$approach1)] <- 0
  allResultsJI$approach2[is.nan(allResultsJI$approach2)] <- 0
  allResultsJI$approach3[is.nan(allResultsJI$approach3)] <- 0
  
  close(pb)
}

#allResultsNDP.ME <- allResultsNDP
allResultsC.ME <-  allResultsC
#allResultsJI.ME <-  allResultsJI


plot(allResultsNDP.ME$Tanimoto_sim, allResultsNDP.ME$approach3, pch = 19, col = "black", xlim = c(0,1), ylim = c(0,1),
     xlab = "Tonimoto similarity", ylab = "spectral similarity")
abline(0,1)
MSE(allResultsNDPME$Tanimoto_sim, allResultsNDPME$approach3)

#plot(allResultsNDP$Tanimoto_dis, allResultsNDP$normal, xlab = "Tanimoto similarity", ylab = "Spectral similarity")

#Boxplot normalized dotproduct
par(fig=c(0,0.33,0.02,1))
boxplot(allResultsNDP$approach1, allResultsNDP$approach2, allResultsNDP$approach3,
        names = c("approach 1", "approach 2", "approach 3"),
        las = 2,
        main = "Normalized dot product",
        ylab = "Similarity index (-)",
        col = c("#5abe7b", "#f8696b", "#638ca6"),
        border = "black",
        ylim = c(0, 1)
)
par(fig=c(0.33,0.66,0.02,1), new=TRUE)
boxplot(allResultsC$approach1, allResultsC$approach2, allResultsC$approach3,
        names = c("approach 1", "approach 2", "approach 3"),
        las = 2,
        main = "Pearson's/Spearman's cor.",
        col = c("#5abe7b", "#f8696b", "#638ca6"),
        border = "black",
        ylim = c(0, 1)
)
par(fig=c(0.66,1,0.02,1), new=TRUE)
boxplot(allResultsNDP$approach1, allResultsNDP$approach2, allResultsNDP$approach3,
        names = c("approach 1", "approach 2", "approach 3"),
        las = 2,
        main = "Jaccard index",
        col = c("#5abe7b", "#f8696b", "#638ca6"),
        border = "black",
        ylim = c(0, 1)
)


################################################################################
#                             Random pair generator                            #
################################################################################

RawData <- Randomize(RawData)

################################################################################
#                                   Clustering                                 #
################################################################################
ColTP2 <- RawData[,2:3]
ColPC2 <- RawData[,12:13]

library("rcdk")

names(ColTP2) <- c("Compound","SMILES")
names(ColPC2) <- c("Compound","SMILES")

Data <- rbind(ColTP2, ColPC2)

Data <- Data[!duplicated(Data$SMILES),]
mols <- parse.smiles(Data$SMILES)

fps <- lapply(mols, get.fingerprint, type='maccs')

fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')

fp.dist <- 1 - fp.sim

row.names(fp.dist) <- Data$Compound

clustering <- hclust(as.dist(fp.dist))

plot(clustering, main='Hierarchical Clustering of CE', xlab = "compounds", ylab = "structural distance (1-similarity)")

################################################################################
#                                     plot graph                               #
################################################################################

plot.new()
par(fig=c(0,0.75,0,1))
plot.window(xlim = xlim, ylim = c(-125, 125))
ticks <- c(-100, -50, 0, 50, 100)
for (i in 1:length(alignment$mz)) lines(rep(alignment$mz[i], 
                                            2), c(0, alignment$intensity.top[i]), col = "#638ca6", lwd = 2)
for (i in 1:length(alignment$mz)) lines(rep(alignment$mz[i], 
                                            2), c(0, -alignment$intensity.bottom[i]), col = "#f8696b", lwd = 2)
axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], 
     ylab = "intensity")
axis(1, pos = -125)
lines(xlim, c(0, 0))
rect(xlim[1], -125, xlim[2], 125)
mtext("m/z", side = 1, line = 2)
mtext("intensity (%)", side = 2, line = 2)
plot.window(xlim = c(0, 20), ylim = c(-10, 10))
text(10, 9, top.label)
text(10, -9, bottom.label)
par(fig=c(0.70,1,0.15,0.95), new=TRUE)
bp <- barplot(height=data$normal, names=data$algorithm, 
              las = 2,
              col="#969696",
              ylab="Similarity (-)", 
              ylim=c(0,1))
text(bp, 0, round(data[,2], 2),cex=0.8,pos=3)

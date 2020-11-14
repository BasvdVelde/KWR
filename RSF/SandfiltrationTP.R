## Script automatically generated on Fri Sep 25 10:33:57 2020

library(patRoon)

#devtools::unload("patRoon")
#remotes::install_github("rickhelmus/patRoon@TP")

# -------------------------
# initialization
# -------------------------


workPath <- "C:/Users/Bas/Google Drive/3. Master/Major stage (KWR)/KWR data/Sand filtration"
#workPath <- "C:\\Users\\veldeba\\Google Drive\\3. Master\\Major stage (KWR)\\KWR data\\Sand filtration"
setwd(workPath)

# Load analysis table
anaInfo <- read.csv("analyses.csv")

# Set to FALSE to skip data pre-treatment
doDataPretreatment <- TRUE
if (doDataPretreatment)
{
  convertMSFiles(anaInfo = anaInfo, from = "thermo",
                 to = "mzML", algorithm = "pwiz", centroid = "vendor")
} 


# -------------------------
# features
# -------------------------


# Find all features.
# NOTE: see manual for many more options
fList <- findFeatures(anaInfo, "openms", noiseThrInt = 4E3,
                      chromFWHM = 3, minFWHM = 1, maxFWHM = 30,
                      chromSNR = 5, mzPPM = 5)

# Group and align features between analysis
fGroups <- groupFeatures(fList, "openms")

# Basic rule based filtering
fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 10000,
                  relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                  blankThreshold = 5, removeBlanks = TRUE,
                  retentionRange = c(138, Inf), mzRange = NULL)

# -----------------------
# Influent and effluent
# ----------------------

FC <- getFCParams(c("influent", "effluent"), thresholdFC = 1)
FC$PVTestFunc <- function(x, y) {
  # print(sprintf("x/y: %f/%f\n", x, y))
  x[x==0] <- runif(sum(x==0), 0, 1)
  y[y==0] <- runif(sum(y==0), 0, 1)
  #x <- log10(x); y <- log10(y)
  
  return(t.test(x, y)$p.value)
}

plotVolcano(fGroups, FC, ylim = c(0,3), averageFunc = median, cex.axis = 1.5, cex.lab = 1.5)


fGroupsMSMS2 <- fGroups[, rGroups = c("influent", "effluent")]
compMSMS <- generateComponentsSpecClust(fGroupsMSMS2, mslists, simMethod = "cosine",
                                        shift = "both", removePrecursor = TRUE, absMzDev = 0.002)
predictTPsComponents(fGroupsMSMS2[, tbl[classification == "decrease"]$group],
                     fGroupsMSMS2[, tbl[classification == "increase"]$group],
                     compMSMS) #geeft error

tbl <- as.data.table(fGroups, FCParams = FC)

tblInc <- tbl[classification == "increase"]

tblInc[tblInc$group == tblInc$group,]

fGroupsInc <- fGroups[, tblInc$group]

# -------------------
# spike and no-spike
# -------------------

FC <- getFCParams(c("effluent_SW", "effluent"), thresholdFC = 2)
FC$PVTestFunc <- function(x, y) {
  # print(sprintf("x/y: %f/%f\n", x, y))
  x[x==0] <- runif(sum(x==0), 0, 1)
  y[y==0] <- runif(sum(y==0), 0, 1)
  #x <- log10(x); y <- log10(y)
  
  return(t.test(x, y, paired = TRUE)$p.value)
}

plotVolcano(fGroups, FC, ylim = c(0,3), averageFunc = median)

tblSW <- as.data.table(fGroups, FCParams = FC)

tblIncSW <- tblSW[classification == "increase"]

tblIncSW[tblIncSW$group == tblIncSW$group,]

fGroupsIncSW <- fGroups[, tblIncSW$group]

metolachlorTPs <- intersect(tblIncSW$group, tblInc$group)
metolachlorTPs <- fGroups[, metolachlorTPs]


susps <- data.table::data.table(name = "metolachlor", SMILES = "CCC1=CC=CC(=C1N(C(C)COC)C(=O)CCl)C")
fGroupsSusp <- groupFeaturesScreening(fGroups, susps, mzWindow = 0.002, adduct = "[M+H]+")

si <- patRoon:::screenInfo(fGroupsSusp)

tbl <- as.data.table(fGroupsSusp, onlyHits = TRUE, collapseSuspects = ",")

fGroupsSuspF <- filter(fGroupsSusp, onlyHits = TRUE)

fGroupsMSMS <- metolachlorTPs[, rGroups = c("effluent_SW", "effluent")]

IDRules <- patRoon:::defaultIDLevelRules()
IDRules[3, "higherThanNext"] <- 0.6
IDRules[IDRules$score == "combScore" | IDRules$score == "isoScore", "relative"] <- FALSE # UNDONE: change in patRoon?

compoundsPos <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsSuspF, mslists, "metfrag", adduct = "[M+H]+",
                                                                                     database = "comptox", maxCandidatesToStop = 1000,
                                                                                     scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", "comptox", onlyDefault = TRUE)$name)))


formulasPos <- withr::with_options(list(patRoon.cache.mode="none"),generateFormulas(fGroupsSuspF, "genform", mslists, relMzDev = 5,
                                                                                    adduct = "[M+H]+", elements = "CHNOPCl",
                                                                                    calculateFeatures = TRUE, featThreshold = 0.75))

fGroupsSuspF <- patRoon:::annotateSuspects(fGroupsSuspF, MSPeakLists = mslists, formulas = formulasPos,
                                               compounds = compoundsPos, IDLevelRules = IDRules,
                                               relMinMSMSIntensity = 0.05)

reportHTML(fGroupsSuspF, compounds = compoundsPos, MSPeakLists = mslists, formulas = formulasPos, path = "SuspPos")

# generate components based on spectral similarity
compMSMS <- generateComponentsSpecClust(fGroupsMSMS, mslists[["M284_R1134_8357"]], simMethod = "cosine", shift = "both", removePrecursor = TRUE,
                                        mzWeight = 0.0, intWeight = 1, absMzDev = 0.002,
                                        relMinIntensity = 0.01)

metolachlorFtrs <- intersect(tblIncSW$group, tblInc$group)

diff <- tblIncSW[tblIncSW$group == "M178_R1440_5263",]$mz-284.1407

patRoon:::specSimilarity(mslists[["M284_R1134_8357"]]$MSMS, mslists[["M178_R1440_5263"]]$MSMS,
                          method = "cosine", shift = "none", precDiff = diff, mzWeight = 0.0, intWeight = 1.0, absMzDev = 0.005, relMinIntensity = 0.1, removePrecursor = TRUE)


#----------------
# Biotransformer
#----------------

TPsBTPos <- withr::with_options(list(patRoon.cache.mode="none"),predictTPsBioTransformer(susps))
TPsBTPos <- filter(TPsBTPos, removeEqualFormulas = TRUE)
#TPs <- filter(TPs, removeEqualFormulas = TRUE, minSimilarity = 0.5)

convertToMFDB(TPsBTPos, "TP-BT.csv", includePrec = TRUE)

fGroupsSuspTPsBTPos <- groupFeaturesScreening(fGroups, convertToSuspects(TPsBTPos), adduct = "[M+H]+", mzWindow = 0.003)
fGroupsSuspTPsBTPos <- filter(fGroupsSuspTPsBTPos, onlyHits=T)
BTFtrs <- as.data.table(fGroupsSuspTPsBTPos)$group

fGroupsSuspTPsBTPosFC <- (as.data.table(fGroupsSuspTPsBTPos, FCParams = FC))
fGroupsSuspTPsBTPosFC[classification == "increase",]

as.data.table(fGroupsSuspTPsBTPos, FCParams = FC)

compoundsTPsBTPos <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsSuspTPsBTPos, mslists, "metfrag", adduct = "[M+H]+",
                                                                                           database = "csv", extraOpts = list(LocalDatabasePath = "TP-BT.csv"),
                                                                                           dbRelMzDev = 10,
                                                                                           fragRelMzDev = 10,
                                                                                           fragAbsMzDev = 0.002,
                                                                                           scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", database = "", onlyDefault = TRUE)$name)))


formulasTPsBTPos <- generateFormulas(fGroupsSuspTPsBTPos, "genform", mslists, relMzDev = 5,
                                     adduct = "[M+H]+", elements = "CHNOPCl",
                                     calculateFeatures = TRUE, featThreshold = 0.75)

fGroupsSuspTPsBTPos <- patRoon:::annotateSuspects(fGroupsSuspTPsBTPos, MSPeakLists = mslists, formulas = formulasTPsBTPos,
                                                  compounds = compoundsTPsBTPos, IDLevelRules = IDRules,
                                                  relMinMSMSIntensity = 0.05)



reportHTML(fGroupsSuspTPsBTPos, compounds = compoundsTPsBTPos, MSPeakLists = mslists, formulas = formulasTPsBTPos, path = "reportRSF")






fGroupsSuspTP <- fGroupsSuspTP[, groupNames(mslistsTP)]
BTFtrs <- as.data.table(fGroupsSuspTP)$group

tblTP <- as.data.table(fGroupsSuspTP, FCParams = FC)
compsTPBT <- generateComponentsTPs(fGroupsSuspTP, fGroupsSuspTP[, tblTP[classification == "increase"]$group], TPs,
                                 mslists, minRTDiff = 0, simMethod = "cosine", removePrecursor = TRUE, mzWeight = 0, 
                                 intWeight = 1, absMzDev = 0.005,  relMinIntensity = 0.1)

compsTPBT[[1]]
componentInfo(compsTP)

plotSpec(compoundsTPsBT, index = 1, groupName = "M178_R1440_5263", MSPeakLists = mslistsTP, xlim = c(50,240))
# annotatedPeakList(compoundsTPsBT, index = 1, groupName = "M250_R1029_7365", MSPeakLists = mslistsTP)
annotatedPeakList(compoundsTPsBT, index = 1, groupName = "M250_R1029_7365", MSPeakLists = mslistsTP, onlyAnnotated=T)
msmsTP <- annotatedPeakList(compoundsTPsBT, index = 1, groupName = "M178_R1440_5263", MSPeakLists = mslistsTP, onlyAnnotated=T)[,1:3]

patRoon:::specSimilarity(msmsTP, mslistsTP[["M178_R1440_5263"]]$MSMS,
                         method = "cosine", shift = "none", precDiff = 0, mzWeight = 0.4, intWeight = 0.1, absMzDev = 0.005, relMinIntensity = 0.01, removePrecursor = FALSE)
#----------------
# metabolic logic
#----------------

TPsML <- predictTPsLogic(fGroupsSuspF, adduct = "[M+H]+", minMass = 40)


fGroupsSuspTPML <- groupFeaturesScreening(fGroups, convertToSuspects(TPsLog), adduct = "[M+H]+", mzWindow = 0.002)
fGroupsSuspTPML <- filter(fGroupsSuspTPLog, onlyHits=T)

MLFtrs <- as.data.table(fGroupsSuspTPML)$group

tblTP <- as.data.table(fGroupsSuspTPML, FCParams = FC)
compsTPML <- generateComponentsTPs(fGroupsSuspTPML, fGroupsSuspTPLog[, tblTP[classification == "increase"]$group], TPsLog,
                                 mslists, minRTDiff = 0, simMethod = "cosine", removePrecursor = TRUE, mzWeight = 0, 
                                 intWeight = 1, absMzDev = 0.002,  relMinIntensity = 0.01)

formsTPsML <- generateFormulas(fGroupsSuspTPML, "genform", MSPeakLists = mslists, elements = "CHNOPSClBr")
compsTPML <- compsTPML[, groupNames(formsTPsML)] 
compsTPML <- filter(compsTPML, formulas = formsTPsML)


compsTPML[[1]]
componentInfo(compsTPML)

fGroupsInteresting <- fGroupsSuspTP[, groupNames(compsTPML)] # no parent!

convertToMFDB(TPs, "test.csv", includePrec = T)

compoundsTPs <- generateCompounds(fGroupsInteresting, mslists, "metfrag", adduct = "[M+H]+",
                                  database = "csv", extraOpts = list(LocalDatabasePath = "NORMAN.csv"))

reportHTML(as(fGroupsInteresting, "featureGroups"), compounds = compoundsTPs, MSPeakLists = mslists, components = NULL)

# -------------------------
# annotation
# -------------------------


# Retrieve MS peak lists
avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
mslists <- generateMSPeakLists(fGroups, "mzr", maxMSRtWindow = 5, precursorMzWindow = 0.5,
                               avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)

mslists <- filter(mslists, withMSMS = TRUE, relMSMSIntThr = 0.1)

mslistsTPs <- generateMSPeakLists(fGroupsSuspTP, "mzr", maxMSRtWindow = 5, precursorMzWindow = 0.5,
                               avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)
mslistsTPs <- filter(mslistsTPs, withMSMS=T)

IDRules <- patRoon:::defaultIDLevelRules()
IDRules[3, "higherThanNext"] <- 0.6

IDRules[IDRules$score == "combScore" | IDRules$score == "isoScore", "relative"] <- FALSE # UNDONE: change in patRoon?

compoundsPos <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsSuspTP, mslistsTPs, "metfrag", adduct = "[M+H]+",
                                                                                      database = "csv", extraOpts = list(LocalDatabasePath = "NORMAN.csv",
                                                                                                                         scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", "csv", onlyDefault = TRUE)$name))))


formulasPos <- generateFormulas(fGroupsSuspTP, "sirius", mslistsTPs, relMzDev = 5,
                                adduct = "[M+H]+", elements = "CHNOPCl",
                                profile = "orbitrap", calculateFeatures = TRUE, featThreshold = 0.75)

fGroupsSuspFPos <- patRoon:::annotateSuspects(fGroupsSuspTP, MSPeakLists = mslistsTPs, formulas = formulasPos,
                                              compounds = compoundsPos, IDLevelRules = IDRules,
                                              relMinMSMSIntensity = 0.05)

reportHTML(fGroupsSuspFPos, compounds = compoundsPos, MSPeakLists = mslistsTPs, formulas = formulasPos)

# uncomment and configure for extra filtering of MS peak lists
# mslists <- filter(mslists, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL,
#                  relMSMSIntThr = NULL, topMSPeaks = NULL, topMSMSPeaks = NULL,
#                  deIsotopeMS = FALSE, deIsotopeMSMS = FALSE)


# Calculate formula candidates
formulas <- generateFormulas(metolachlorTPs, "sirius", mslists, relMzDev = 5,
                             adduct = "[M+H]+", elements = "CHNOPCl",
                             profile = "orbitrap", calculateFeatures = TRUE, featThreshold = 0.75)

formulas <- withr::with_options(list(patRoon.cache.mode="none"), generateFormulas(metolachlorTPs, "genform", mslistsTP, relMzDev = 5,
                             adduct = "[M+H]+", elements = "CHNOPCl",
                             calculateFeatures = TRUE, featThreshold = 0.75))

# Find compound structure candidates
compounds <- generateCompounds(metolachlorTPs, mslists, "metfrag", relMzDev = 5,
                               adduct = "[M+H]+", elements = "CHNOPCl", profile = "orbitrap",
                               metfragDatabase = "pubchem")
View(as.data.table(compounds))

compounds <- withr::with_options(list(patRoon.cache.mode="none"), generateCompounds(fGroupsSuspTP, mslists, "metfrag", method = "CL", dbRelMzDev = 5,
                               fragRelMzDev = 5, fragAbsMzDev = 0.002,
                               adduct = "[M+H]+", database = "pubchemlite", maxCandidatesToStop = 300))


# -------------------------
# reporting
# -------------------------


reportCSV(metolachlorTPs, path = "report", reportFeatures = FALSE, formulas = formulas,
          compounds = compounds, compoundsNormalizeScores = "max",
          components = NULL)

reportHTML(metolachlorTPs, path = "report", reportPlots = c("chord", "venn", "upset", "eics", "formulas"),
           formulas = formulas, compounds = compounds, compoundsNormalizeScores = "max",
           components = NULL, MSPeakLists = mslistsTP,
           selfContained = FALSE, openReport = TRUE)


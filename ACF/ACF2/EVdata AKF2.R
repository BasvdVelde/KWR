# --------------------------
#  Positive ionization mode
# --------------------------

library(patRoon)

#devtools::unload("patRoon")
#remotes::install_github("rickhelmus/patRoon@TP")
#remotes::install_github("rickhelmus/patRoon")

# -------------------------
# initialization Pos
# -------------------------

workPath <- "C:/Users/Bas/Google Drive/3. Master/Major stage (KWR)/KWR data/EV data AKF2"
setwd(workPath)

# Load analysis table
anaInfoPos <- read.csv("analysesPosAKF2.csv")

# Set to FALSE to skip data pre-treatment
doDataPretreatment <- TRUE
if (doDataPretreatment)
{
  convertMSFiles(anaInfo = anaInfoPos, from = "thermo",
                 to = "mzML", algorithm = "pwiz", centroid = "vendor")
} 

# -------------------------
# features Pos
# -------------------------

fListPos <- findFeatures(anaInfoPos, "openms", noiseThrInt = 4E3,
                         chromFWHM = 3, minFWHM = 1, maxFWHM = 30,
                         chromSNR = 5, mzPPM = 20)

fGroupsPos <- groupFeatures(fListPos, "openms")

fGroupsPos <- filter(fGroupsPos, preAbsMinIntensity = 100, absMinIntensity = 10000,
                     relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                     blankThreshold = 5, removeBlanks = TRUE,
                     retentionRange = c(120, Inf), mzRange = NULL)

avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)

mslistsPos <- generateMSPeakLists(fGroupsPos, "mzr", maxMSRtWindow = 5, precursorMzWindow = 0.5,
                                  avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)

mslistsPos <- filter(mslistsPos, withMSMS = TRUE, relMSMSIntThr = 0.05)

# ----------------------
# Suspect screening Pos
# ----------------------

FC <- getFCParams(c("influent-pos", "effluent-pos"), thresholdFC = 1)
FC$PVTestFunc <- function(x, y) {
  # print(sprintf("x/y: %f/%f\n", x, y))
  x[x==0] <- runif(sum(x==0), 0, 1)
  y[y==0] <- runif(sum(y==0), 0, 1)
  #x <- log10(x); y <- log10(y)
  
  return(t.test(x, y, paired = TRUE)$p.value)
}

plotVolcano(fGroupsPos, FC, ylim = c(0,3), averageFunc = median)

tblPos <- as.data.table(fGroupsPos, FCParams = FC)

tblIncPos <- tblPos[classification == "increase"]
write.csv(tblIncPos, "Inc.csv")
tblDecPos <- tblPos[classification == "decrease"]

IncFtrs <- tblIncPos$group

tblIncPos[tblIncPos$group == tblIncPos$group,]

fGroupsIncPos <- fGroupsPos[, tblIncPos$group]

TPsPos <- intersect(tblIncPos$group, tblPos$group)


susps <- read.csv("EVsusp.csv")

fGroupsSuspPos <- groupFeaturesScreening(fGroupsPos, susps, mzWindow = 0.005, adduct = "[M+H]+")

siPos <- patRoon:::screenInfo(fGroupsSuspPos)

tblPos <- as.data.table(fGroupsSuspPos, onlyHits = TRUE, collapseSuspects = ",")
tblPos1 <-  tblPos[,c(1:3,10)]
suspGroupPos <- tblPos1$group

fGroupsSuspFPos <- filter(fGroupsSuspPos, onlyHits = TRUE)

IDRules <- patRoon:::defaultIDLevelRules()
IDRules[3, "higherThanNext"] <- 0.6
IDRules[IDRules$score == "combScore" | IDRules$score == "isoScore", "relative"] <- FALSE # UNDONE: change in patRoon?



fGroupsMSMS <- TPsPos[, rGroups = c("effluent-pos", "effluent-pos")]


compoundsPos <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsSuspFPos, mslistsPos, "metfrag", adduct = "[M+H]+",
                                                                                      database = "comptox", maxCandidatesToStop = 1000,
                                                                                      scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", "comptox", onlyDefault = TRUE)$name)))

compoundsIncPos <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsIncPos, mslistsPos, "metfrag", adduct = "[M+H]+",
                                                                                      database = "comptox", maxCandidatesToStop = 3000,
                                                                                      scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", "comptox", onlyDefault = TRUE)$name)))


formulasPos <- withr::with_options(list(patRoon.cache.mode="none"),generateFormulas(fGroupsSuspFPos, "sirius", mslistsPos, relMzDev = 5,
                                adduct = "[M+H]+", elements = "CHNOPCl",
                                profile = "orbitrap", calculateFeatures = TRUE, featThreshold = 0.75))

formulasIncPos <- generateFormulas(fGroupsIncPos, "genform", mslistsPos, relMzDev = 5,
                                adduct = "[M+H]+", elements = "CHNOCl",
                                calculateFeatures = TRUE, featThreshold = 0.75)

fGroupsSuspFPos1 <- patRoon:::annotateSuspects(fGroupsSuspFPos, MSPeakLists = mslistsPos, formulas = formulasPos,
                                              compounds = compoundsPos, IDLevelRules = IDRules,
                                              relMinMSMSIntensity = 0.05)

# fGroupsSuspFPos <- filter(fGroupsSuspFPos, maxLevel = 3, maxCompRank = 10, maxFormRank = 10,
#                           selectBy = "level", onlyHits = TRUE) # selectHitsBy soon



reportHTML(fGroupsSuspFPos, compounds = compoundsPos, MSPeakLists = mslistsPos, formulas = formulasPos, path = "SuspPos")

reportHTML(fGroupsIncPos, compounds = compoundsIncPos, MSPeakLists = mslistsPos, formulas = formulasIncPos, path = "IncPos")
# ----------------------
# Biotransformer Pos
# ----------------------

TPsBTPos <- withr::with_options(list(patRoon.cache.mode="none"),predictTPsBioTransformer(susps))
TPsBTPos <- filter(TPsBTPos, removeEqualFormulas = TRUE)
#TPs <- filter(TPs, removeEqualFormulas = TRUE, minSimilarity = 0.5)

convertToMFDB(TPsBTPos, "TP-BT.csv", includePrec = TRUE)

fGroupsSuspTPsBTPos <- groupFeaturesScreening(fGroupsPos, convertToSuspects(TPsBTPos), adduct = "[M+H]+", mzWindow = 0.003)
fGroupsSuspTPsBTPos <- filter(fGroupsSuspTPsBTPos, onlyHits=T)
BTFtrs <- as.data.table(fGroupsSuspTPsBTPos)$group

fGroupsSuspTPsBTPosFC <- (as.data.table(fGroupsSuspTPsBTPos, FCParams = FC))
fGroupsSuspTPsBTPosFC[classification == "increase",]

as.data.table(fGroupsSuspTPsBTPos, FCParams = FC)

compoundsTPsBTPos <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsSuspTPsBTPos, mslistsPos, "metfrag", adduct = "[M+H]+",
                                                                                           database = "csv", extraOpts = list(LocalDatabasePath = "TP-BT.csv"),
                                                                                           dbRelMzDev = 10,
                                                                                           fragRelMzDev = 10,
                                                                                           fragAbsMzDev = 0.002,
                                                                                           scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", database = "", onlyDefault = TRUE)$name)))


formulasTPsBTPos <- generateFormulas(fGroupsSuspTPsBTPos, "sirius", mslistsPos, relMzDev = 5,
                                     adduct = "[M+H]+", elements = "CHNOPCl",
                                     profile = "orbitrap", calculateFeatures = TRUE, featThreshold = 0.75)

fGroupsSuspTPsBTPos <- patRoon:::annotateSuspects(fGroupsSuspTPsBTPos, MSPeakLists = mslistsPos, formulas = formulasTPsBTPos,
                                                  compounds = compoundsTPsBTPos, IDLevelRules = IDRules,
                                                  relMinMSMSIntensity = 0.05)



reportHTML(fGroupsSuspTPsBTPos, compounds = compoundsTPsBTPos, MSPeakLists = mslistsPos, formulas = formulasTPsBTPos, components = NULL, path = "reportBTAKF2")

#----------------------
# metabolic logic Pos
#----------------------

TPsMLPos <- predictTPsLogic(fGroupsSuspFPos, adduct = "[M+H]+", minMass = 40)


fGroupsSuspTPsMLPos <- groupFeaturesScreening(fGroupsPos, convertToSuspects(TPsMLPos), adduct = "[M+H]+", mzWindow = 0.002)
fGroupsSuspTPsMLPos <- filter(fGroupsSuspTPsMLPos, onlyHits=T)

MLFtrs <- as.data.table(fGroupsSuspTPsMLPos)$group

tblTP <- as.data.table(fGroupsSuspTPsMLPos, FCParams = FC)
compsTPML <- generateComponentsTPs(fGroupsSuspTPML, fGroupsSuspTPML[, tblTP[classification == "increase"]$group], TPsML,
                                   mslistsPos, minRTDiff = 0, simMethod = "cosine", removePrecursor = TRUE, mzWeight = 0, 
                                   intWeight = 1, absMzDev = 0.002,  relMinIntensity = 0.01)

formsTPsML <- generateFormulas(fGroupsSuspTPML, "genform", MSPeakLists = mslistsPos, elements = "CHNOPSClBr")
compsTPML <- compsTPML[, groupNames(formsTPsML)] 
compsTPML <- filter(compsTPML, formulas = formsTPsML)

convertToMFDB(TPsMLPos, "TP-ML.csv", includePrec = TRUE)

compsTPML[[1]]
componentInfo(compsTPML)

fGroupsInteresting <- fGroupsSuspTP[, groupNames(compsTPML)] # no parent!

convertToMFDB(TPs, "test.csv", includePrec = T)


compoundsTPsMLPos <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsSuspTPsMLPos, mslistsPos, "metfrag", adduct = "[M+H]+",
                                                                                           database = "csv", extraOpts = list(LocalDatabasePath = "TP-BT.csv"),
                                                                                           dbRelMzDev = 10,
                                                                                           fragRelMzDev = 10,
                                                                                           fragAbsMzDev = 0.002,
                                                                                           scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", database = "", onlyDefault = TRUE)$name)))


formulasTPsMLPos <- generateFormulas(fGroupsSuspTPsBTPos, "sirius", mslistsPos, relMzDev = 5,
                                     adduct = "[M+H]+", elements = "CHNOPCl",
                                     profile = "orbitrap", calculateFeatures = TRUE, featThreshold = 0.75)

fGroupsSuspTPsMLPos <- patRoon:::annotateSuspects(fGroupsSuspTPsBTPos, MSPeakLists = mslistsPos, formulas = formulasTPsBTPos,
                                                  compounds = compoundsTPsBTPos, IDLevelRules = IDRules,
                                                  relMinMSMSIntensity = 0.05)



reportHTML(fGroupsSuspTPsMLPos, compounds = compoundsTPsBTPos, MSPeakLists = mslistsPos, formulas = formulasTPsBTPos, components = NULL, path = "reportBTAKF2")


diff = (tblPos[tblPos$group == "M166_R135_4777"]$mz - tblPos[tblPos$group == "M182_R135_5609"]$mz)
  
patRoon:::specSimilarity(mslistsPos[["M166_R135_4777"]]$MSMS, mslistsPos[["M182_R135_5609"]]$MSMS,
                           method = "cosine", shift = "both", precDiff = diff, mzWeight = 0.0, intWeight = 1.0, absMzDev = 0.005, relMinIntensity = 0.01, removePrecursor = TRUE)

#----------------------
# NTS Pos
#----------------------

fGroupsMSMSPos <- fGroupsPos[, rGroups = c("influent-pos", "effluent-pos")]

#fGroupsMSMSPos <- fGroupsMSMSPos[, setdiff(names(fGroupsMSMSPos),
#                                           c("M214_R141_6857"))]

# get log2fc data
tblFCPos <- as.data.table(fGroupsMSMSPos, FCParams = getFCParams(replicateGroups(fGroupsMSMSPos)))
fGroupsMSMSPos <- fGroupsMSMSPos[, c(suspGroupPos, tblFCPos[classification == "increase"]$group)]

# generate components based on spectral similarity
compMSMSPos <- generateComponentsSpecClust(fGroupsMSMSPos, mslistsPos, simMethod = "cosine",
                                           shift = "both", removePrecursor = TRUE, absMzDev = 0.005)

# the clusters are automatically assigned, but you most likely want to re-do it manually:
compMSMSPos <- treeCut(compMSMSPos, h = 0.99) # min 0.5 similarity


# make TP predictions based on the simularity clustering and log2fc
TPsMSMSPos <- predictTPsComponents(compMSMSPos, fGroupsMSMSPos[, suspGroupPos],
                                   fGroupsMSMSPos[, tblFCPos[classification == "increase"]$group])

# and finally generate TP components 'as usual'. It probably makes sense to  use
# the same spectral similarity parameters as for generateComponentsSpecClust() :-)
componTPMSMSPos <- generateComponentsTPs(fGroupsMSMSPos, pred = TPsMSMSPos, MSPeakLists = mslistsPos,
                                         simMethod = "cosine", removePrecursor = TRUE,
                                         mzWeight = 0, intWeight = 1, absMzDev = 0.005,
                                         relMinIntensity = 0.05)

tblCompon <- as.data.table(unlist(componTPMSMSPos))

tblCompon <- tblCompon[tblCompon$specSimilarityBoth >= 0.5,]

write.csv(unlist(tblCompon), "tblComponPos.csv")

tblCompon[tblCompon$precursor_group == suspGroupPos]

componTPMSMSPosSub <- componTPMSMSPos[tbl$name, tbl$group]

# preview results
reportHTML(fGroupsMSMSPos, components = componTPMSMSPos, path = "SpecSim")

fGroupsMSMSPos1 <- fGroupsMSMSPos[, union(groupNames(compoundsIncPos), groupNames(formulasIncPos))]

Spec.simPos <- calcMSMSSims(fGroupsMSMSPos[, suspGroupPos], fGroupsMSMSPos1[, tblFCPos[classification == "increase"]$group], mslistsPos)
write.csv(Spec.simPos, "specsimPos.csv")

# --------------------------
#  Negative ionization mode
# --------------------------

#devtools::unload("patRoon")
#remotes::install_github("rickhelmus/patRoon@TP")
#remotes::install_github("rickhelmus/patRoon")

# -------------------------
# initialization Neg
# -------------------------

# Load analysis table
anaInfoNeg <- read.csv("analysesNegAKF2.csv")

# Set to FALSE to skip data pre-treatment
doDataPretreatment <- TRUE
if (doDataPretreatment)
{
  convertMSFiles(anaInfo = anaInfoNeg, from = "thermo",
                 to = "mzML", algorithm = "pwiz", centroid = "vendor")
} 

# -------------------------
# features Neg
# -------------------------

fListNeg <- findFeatures(anaInfoNeg, "openms", noiseThrInt = 4E3,
                         chromFWHM = 3, minFWHM = 1, maxFWHM = 30,
                         chromSNR = 5, mzPPM = 5)

fGroupsNeg <- groupFeatures(fListNeg, "openms")

fGroupsNeg <- filter(fGroupsNeg, preAbsMinIntensity = 100, absMinIntensity = 10000,
                     relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                     blankThreshold = 5, removeBlanks = TRUE,
                     retentionRange = c(120, Inf), mzRange = NULL)

avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)

mslistsNeg <- generateMSPeakLists(fGroupsNeg, "mzr", maxMSRtWindow = 5, precursorMzWindow = 0.5,
                                  avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)

mslistsNeg <- filter(mslistsNeg, withMSMS = TRUE, relMSMSIntThr = 0.05) # extra filter: topMost = 10

# ----------------------
# Suspect screening Neg
# ----------------------

FCNeg <- getFCParams(c("influent-neg", "effluent-neg"), thresholdFC = 1)
FC$PVTestFunc <- function(x, y) {
  # print(sprintf("x/y: %f/%f\n", x, y))
  x[x==0] <- runif(sum(x==0), 0, 1)
  y[y==0] <- runif(sum(y==0), 0, 1)
  #x <- log10(x); y <- log10(y)
  
  return(t.test(x, y, paired = TRUE)$p.value)
}

plotVolcano(fGroupsNeg, FCNeg, ylim = c(0,3), averageFunc = median)

tblNeg <- as.data.table(fGroupsNeg, FCParams = FCNeg)

tblIncNeg <- tblNeg[classification == "increase"]
write.csv(tblIncNeg, "IncNeg.csv")
tblDecNeg <- tblNeg[classification == "decrease"]

IncFtrs <- tblIncNeg$group

tblIncNeg[tblIncNeg$group == tblIncNeg$group,]

fGroupsIncNeg <- fGroupsNeg[, tblIncNeg$group]

TPsNeg <- intersect(tblIncNeg$group, tblNeg$group)


susps <- read.csv("EVsusp.csv")

fGroupsSuspNeg <- groupFeaturesScreening(fGroupsNeg, susps, mzWindow = 0.005, adduct = "[M-H]-")

siNeg <- patRoon:::screenInfo(fGroupsSuspNeg)

tblNeg <- as.data.table(fGroupsSuspNeg, onlyHits = TRUE, collapseSuspects = ",")
tblNeg1 <-  tblNeg[,c(1:3,10)]
suspGroupNeg <- tblNeg1$group

fGroupsSuspFNeg <- filter(fGroupsSuspNeg, onlyHits = TRUE)

IDRules <- patRoon:::defaultIDLevelRules()
IDRules[3, "higherThanNext"] <- 0.6
IDRules[IDRules$score == "combScore" | IDRules$score == "isoScore", "relative"] <- FALSE # UNDONE: change in patRoon?



fGroupsMSMSNeg <- TPsNeg[, rGroups = c("effluent-neg", "effluent-neg")]


compoundsNeg <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsSuspFNeg, mslistsNeg, "metfrag", adduct = "[M-H]-",
                                                                                      database = "comptox", maxCandidatesToStop = 1000,
                                                                                      scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", "comptox", onlyDefault = TRUE)$name)))

compoundsIncNeg <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsIncNeg, mslistsNeg, "metfrag", adduct = "[M-H]-",
                                                                                         database = "comptox", maxCandidatesToStop = 3000,
                                                                                         scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", "pubchem", onlyDefault = TRUE)$name)))


formulasNeg <- generateFormulas(fGroupsSuspFNeg, "genform", mslistsNeg, relMzDev = 5,
                                adduct = "[M-H]-", elements = "CHNOPCl",
                                profile = "orbitrap", calculateFeatures = TRUE, featThreshold = 0.75)

formulasIncNeg <- generateFormulas(fGroupsIncNeg, "genform", mslistsNeg, relMzDev = 5,
                                   adduct = "[M-H]-", elements = "CHNOPCl",
                                   calculateFeatures = TRUE, featThreshold = 0.75)

fGroupsSuspFNeg <- patRoon:::annotateSuspects(fGroupsSuspFNeg, MSPeakLists = mslistsNeg, formulas = formulasNeg,
                                              compounds = compoundsNeg, IDLevelRules = IDRules,
                                              relMinMSMSIntensity = 0.05)

# fGroupsSuspFNeg <- filter(fGroupsSuspFNeg, maxLevel = 3, maxCompRank = 10, maxFormRank = 10,
#                           selectBy = "level", onlyHits = TRUE) # selectHitsBy soon



reportHTML(fGroupsSuspFNeg, compounds = compoundsNeg, MSPeakLists = mslistsNeg, formulas = formulasNeg, path = "SuspNeg")

reportHTML(fGroupsIncNeg, compounds = compoundsIncNeg, MSPeakLists = mslistsNeg, formulas = formulasIncNeg, path = "IncNeg")
# ----------------------
# Biotransformer Neg
# ----------------------

TPsBTNeg <- withr::with_options(list(patRoon.cache.mode="none"),predictTPsBioTransformer(susps))
TPsBTNeg <- filter(TPsBTNeg, removeEqualFormulas = TRUE)
#TPs <- filter(TPs, removeEqualFormulas = TRUE, minSimilarity = 0.5)

convertToMFDB(TPsBTNeg, "TP-BTNeg.csv", includePrec = TRUE)

fGroupsSuspTPsBTNeg <- groupFeaturesScreening(fGroupsNeg, convertToSuspects(TPsBTNeg), adduct = "[M-H]-", mzWindow = 0.005)
fGroupsSuspTPsBTNeg <- filter(fGroupsSuspTPsBTNeg, onlyHits=T)
BTFtrsNeg <- as.data.table(fGroupsSuspTPsBTNeg)$group

fGroupsSuspTPsBTNegFC <- (as.data.table(fGroupsSuspTPsBTNeg, FCParams = FC))
fGroupsSuspTPsBTNegFC[classification == "increase",]

as.data.table(fGroupsSuspTPsBTNeg, FCParams = FC)

compoundsTPsBTNeg <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsSuspTPsBTNeg, mslistsNeg, "metfrag", adduct = "[M-H]-",
                                                                                           database = "csv", extraOpts = list(LocalDatabasePath = "TP-BTNeg.csv"),
                                                                                           dbRelMzDev = 10,
                                                                                           fragRelMzDev = 10,
                                                                                           fragAbsMzDev = 0.002,
                                                                                           scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", database = "", onlyDefault = TRUE)$name)))


formulasTPsBTNeg <- generateFormulas(fGroupsSuspTPsBTNeg, "sirius", mslistsNeg, relMzDev = 5,
                                     adduct = "[M-H]-", elements = "CHNOPCl",
                                     profile = "orbitrap", calculateFeatures = TRUE, featThreshold = 0.75)

fGroupsSuspTPsBTNeg <- patRoon:::annotateSuspects(fGroupsSuspTPsBTNeg, MSPeakLists = mslistsNeg, formulas = formulasTPsBTNeg,
                                                  compounds = compoundsTPsBTNeg, IDLevelRules = IDRules,
                                                  relMinMSMSIntensity = 0.05)



reportHTML(fGroupsSuspTPsBTNeg, compounds = compoundsTPsBTNeg, MSPeakLists = mslistsNeg, formulas = formulasTPsBTNeg, components = NULL, path = "reportBTAKF2Neg")

#----------------------
# metabolic logic Neg
#----------------------

TPsMLNeg <- predictTPsLogic(fGroupsSuspFNeg, adduct = "[M+H]+", minMass = 40)


fGroupsSuspTPsMLNeg <- groupFeaturesScreening(fGroupsNeg, convertToSuspects(TPsMLNeg), adduct = "[M+H]+", mzWindow = 0.002)
fGroupsSuspTPsMLNeg <- filter(fGroupsSuspTPsMLNeg, onlyHits=T)

MLFtrsNeg <- as.data.table(fGroupsSuspTPsMLNeg)$group

tblTP <- as.data.table(fGroupsSuspTPsMLNeg, FCParams = FC)
compsTPML <- generateComponentsTPs(fGroupsSuspTPML, fGroupsSuspTPML[, tblTP[classification == "increase"]$group], TPsML,
                                   mslistsNeg, minRTDiff = 0, simMethod = "cosine", removePrecursor = TRUE, mzWeight = 0, 
                                   intWeight = 1, absMzDev = 0.002,  relMinIntensity = 0.01)

formsTPsML <- generateFormulas(fGroupsSuspTPML, "genform", MSPeakLists = mslistsNeg, elements = "CHNOPSClBr")
compsTPML <- compsTPML[, groupNames(formsTPsML)] 
compsTPML <- filter(compsTPML, formulas = formsTPsML)

convertToMFDB(TPsMLNeg, "TP-ML.csv", includePrec = TRUE)

compsTPML[[1]]
componentInfo(compsTPML)

fGroupsInteresting <- fGroupsSuspTP[, groupNames(compsTPML)] # no parent!

convertToMFDB(TPs, "test.csv", includePrec = T)


compoundsTPsMLNeg <- withr::with_options(list(patRoon.cache.mode="none"),generateCompounds(fGroupsSuspTPsMLNeg, mslistsNeg, "metfrag", adduct = "[M+H]+",
                                                                                           database = "csv", extraOpts = list(LocalDatabasePath = "TP-BT.csv"),
                                                                                           dbRelMzDev = 10,
                                                                                           fragRelMzDev = 10,
                                                                                           fragAbsMzDev = 0.002,
                                                                                           scoreTypes = c("individualMoNAScore", compoundScorings("metfrag", database = "", onlyDefault = TRUE)$name)))


formulasTPsMLNeg <- generateFormulas(fGroupsSuspTPsBTNeg, "sirius", mslistsNeg, relMzDev = 5,
                                     adduct = "[M+H]+", elements = "CHNOPCl",
                                     profile = "orbitrap", calculateFeatures = TRUE, featThreshold = 0.75)

fGroupsSuspTPsMLNeg <- patRoon:::annotateSuspects(fGroupsSuspTPsBTNeg, MSPeakLists = mslistsNeg, formulas = formulasTPsBTNeg,
                                                  compounds = compoundsTPsBTNeg, IDLevelRules = IDRules,
                                                  relMinMSMSIntensity = 0.05)



reportHTML(fGroupsSuspTPsMLNeg, compounds = compoundsTPsBTNeg, MSPeakLists = mslistsNeg, formulas = formulasTPsBTNeg, components = NULL, path = "reportBTAKF2")


diff = (tblNeg[tblNeg$group == "M285_R851_3826"]$mz - tblNeg[tblNeg$group == "M182_R135_5609"]$mz)

patRoon:::specSimilarity(mslistsNeg[["M285_R851_3826"]]$MSMS, mslistsNeg[["M182_R135_5609"]]$MSMS,
                         method = "cosine", shift = "both", precDiff = diff, mzWeight = 0.0, intWeight = 1.0, absMzDev = 0.005, relMinIntensity = 0.01, removePrecursor = TRUE)


#----------------------
# NTS Neg
#----------------------

fGroupsMSMSNeg <- fGroupsNeg[, rGroups = c("influent-neg", "effluent-neg")]

# generate components based on spectral similarity
compMSMSNeg <- generateComponentsSpecClust(fGroupsMSMSNeg, mslistsNeg, simMethod = "cosine",
                                           shift = "both", removePrecursor = TRUE, absMzDev = 0.002)

# the clusters are automatically assigned, but you most likely want to re-do it manually:
compMSMSNeg <- treeCut(compMSMSNeg, h = 0.5) # min 0.5 similarity

# get log2fc data
tblFCNeg <- as.data.table(fGroupsMSMSNeg, FCParams = getFCParams(replicateGroups(fGroupsMSMSNeg)))

# make TP predictions based on the simularity clustering and log2fc
TPsMSMSNeg <- predictTPsComponents(fGroupsMSMSNeg[, tblFCNeg[classification == "decrease"]$group],
                                   fGroupsMSMSNeg[, tblFCNeg[classification == "increase"]$group],
                                   compMSMSNeg)

# and finally generate TP components 'as usual'. It probably makes sense to  use
# the same spectral similarity parameters as for generateComponentsSpecClust() :-)
componTPMSMSNeg <- generateComponentsTPs(fGroupsMSMSNeg, pred = TPsMSMSNeg, MSPeakLists = mslistsNeg,
                                         simMethod = "cosine", removePrecursor = TRUE,
                                         mzWeight = 0, intWeight = 1, absMzDev = 0.002,
                                         relMinIntensity = 0.1)

# preview results
reportHTML(fGroupsMSMS, components = componTPMSMS)

fGroupsMSMSNeg1 <- fGroupsMSMSNeg[, union(groupNames(compoundsIncNeg), groupNames(formulasIncNeg))]

Spec.simNeg <- calcMSMSSims(fGroupsMSMSNeg[, suspGroupNeg], fGroupsMSMSNeg1[, tblFCNeg[classification == "increase"]$group], mslistsNeg)
write.csv(Spec.simNeg, "specsimNeg.csv")

calcMSMSSims <- function(fGroupsParent, fGroupsTPs, mslists)
{
  ret <- data.table::rbindlist(lapply(names(fGroupsTPs), function(grpTP)
  {
    if (is.null(mslists[[grpTP]]) || is.null(mslists[[grpTP]][["MSMS"]]))
      return(data.table::data.table())
    sapply(names(fGroupsParent), function(grpPar)
    {
      if (is.null(mslists[[grpPar]]) || is.null(mslists[[grpPar]][["MSMS"]]))
        return(NA)
      precDiff <- groupInfo(fGroupsParent)[grpPar, "mzs"] - groupInfo(fGroupsTPs)[grpTP, "mzs"]
      return(patRoon:::specSimilarity(mslists[[grpPar]]$MSMS, mslists[[grpTP]]$MSMS,
                                      "cosine", shift = "both", precDiff = precDiff, removePrecursor = TRUE,
                                      mzWeight = 0.4, intWeight = 0.1, absMzDev = 0.005, relMinIntensity = 0.05))
    }, simplify = FALSE)
  }))
  ret <- as.data.frame(ret); rownames(ret) <- names(fGroupsTPs)
  return(ret)
}

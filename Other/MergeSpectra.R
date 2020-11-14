library("rjson")
library("data.table")
source("utils.R")
RawData <- read.csv(url("https://drive.google.com/uc?export=download&id=1m-_JA_XDw1LYofSyyHg2br4Kmvo1OtYm")) 

for (i in seq_len(nrow(RawData))) {

DF <- apply(RawData[i, 5:10], 2, function(x) MBID(x))
precursor <- apply(RawData[i, 5:10], 2, function(x) MBID(x, precursor = TRUE))
spcomb <- as.data.table(do.call(rbind.data.frame, DF))
h=0.005
b=0.0
{ 
setorder(spcomb, mz)

mzd <- dist(spcomb$mz)
hc <- fastcluster::hclust(mzd)
spcomb[, cluster := cutree(hc, h = h)]

spcount <- length(spcomb)
ret <- spcomb[, .(mz = mean(mz), intensity = sum(intensity) / spcount), by = cluster][, cluster := NULL][]

ret <- ret[intensity >= b]

ret[, intensity := intensity / max(intensity) * 100]

} 

peaks <- ret
n_peaks <- nrow(ret)
precursor <- unlist(precursor)
precursor_mz <- max(precursor)
lst <- list(n_peaks, peaks, precursor_mz)
names(lst) <- c("n_peaks", "peaks", "precursor_mz")

spl <- split(lst$peaks, f = seq_len(nrow(lst$peaks)))
names(spl) <- NULL

spl <- lapply(spl, unlist, use.names = FALSE)

lst <- list(n_peaks, spl, precursor_mz)
names(lst) <- c("n_peaks", "peaks", "precursor_mz")
jsonData <- toJSON(lst)
write(jsonData, paste0("Merge2/METP", RawData$TP_ID[i], ".json"))

}



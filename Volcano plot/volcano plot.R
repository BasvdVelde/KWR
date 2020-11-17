#collect the data into a matrix
b_a <- as.data.frame(fGroups, average = TRUE) #generates the average data for the two data sets

group <- as.data.frame(b_a["group"]) #generated the features in the MS data

datap <- as.data.frame(fGroups)


datapp <- datap[4:9]

datappp <- as.data.frame(datap[, c(1, 4:9)])

for(i in seq_along(datapp)) datapp[[i]][datapp[[i]] == 0] <- runif(length(i), 0, 1)

before <- colnames(b_a[4]) 

after <- colnames(b_a[5])

before_data <- apply(as.data.frame(datapp[,1:3]), 1, FUN = median)
after_data <- apply(as.data.frame(datapp[,4:6]), 1, FUN = median)

FC <- (after_data/before_data)
logFC <- log2(FC)

pv <- sapply(1:nrow(datapp), function(i) t.test(as.numeric(datapp[i, 1:3]), as.numeric(datapp[i, 4:6]))$p.value)

pvg <- cbind(group, pv)

pvc <- p.adjust(pv, method = "BH", n = length(pv))

logPV <- -log10(pvc)

resultTable <- data.frame(group = pvg$group, FC = logFC, FCSW = logFCSW, pvc = logPV, color = "#969696")
isSignificant <- resultTable$pvc > 1.30103
isTP <- isSignificant == TRUE & resultTable$FC > 1
isTPfromSpikeIn <- isSignificant == TRUE & resultTable$FC > 1 & resultTable$FCSW > 2
isPC <- isSignificant & resultTable$FC < 1
signFC <- !isSignificant & abs(resultTable$FC) > 1
signPV <- isSignificant & abs(resultTable$FC) < 1

resultTable$color[isTP] <- "#5abe7b"
resultTable$color[isTPfromSpikeIn] <- "#3d8555"
resultTable$color[isPC] <- "#f8696b"
resultTable$color[signFC] <- "#638ca6"
resultTable$color[signPV] <- "#fff6dd"

#log2FC >2 extra color
datappSW <- datap[7:12]
for(i in seq_along(datappSW)) datappSW[[i]][datappSW[[i]] == 0] <- runif(length(i), 0, 1)

before_dataSW <- apply(as.data.frame(datappSW[,1:3]), 1, FUN = median)
after_dataSW <- apply(as.data.frame(datappSW[,4:6]), 1, FUN = median)

FCSW <- (before_dataSW/after_dataSW)
logFCSW <- log2(FCSW)


#Make volcano plot
plot(resultTable$FC, resultTable$pv, ylab = "-log10 p value", 
     xlab = "log2 fold change", col = resultTable$color, pch = 19, ylim = c(0,4))
#text(resultTable$FC, resultTable$pv, labels = resultTable$group , pos = 4, cex = 0.6)
abline(v = c(-1, 1), col = c("red", "red"), lty = c(2, 2), lwd = c(2, 2), h = -log10(0.05))
legend("topright", legend = c("Possible TP", "Possible PC", "Significant FC", "Significant PV", "Not significant", "TP from spikein"),
       col=c("#5abe7b", "#f8696b", "#638ca6", "#fff6dd", "#969696", "#3d8555"), cex = 0.8, pch = 19)

#collect MS/MS data for the TPs and PCs
groupsTP <- resultTable$group[isTP]

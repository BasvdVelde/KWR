RawData <- read.csv(url("https://drive.google.com/uc?export=download&id=1m-_JA_XDw1LYofSyyHg2br4Kmvo1OtYm")) 
url <- "https://metabolomics-usi.ucsd.edu/json/?usi=mzspec:MASSBANK:"

PC <- RawData[,5:10]

names(PC) <- c("CE15", "CE30", "CE45", "CE60", "CE75" , "CE90")

TP <- RawData[,15:20]

names(TP) <- c("CE15", "CE30", "CE45", "CE60", "CE75" , "CE90")

RawData1 <- rbind(PC, TP)

filenames <- c(RawData1$CE15, RawData1$CE30, RawData1$CE45, RawData1$CE60, RawData1$CE75, RawData1$CE90)
filenames <- filenames[!is.na(filenames)]


for (i in 1:length(filenames)) {
  download.file(paste(url, filenames[i], sep = ""), paste(getwd(), "/", filenames[i], ".JSON",
                                                      sep = ""))
}

setwd("C:\\Users\\Bas\\Desktop\\MBID")


filelist <- list.files(pattern = ".json")
someVar <- lapply(filelist, function(x) { 
  textfile <- fromJSON(x)
  write.csv(textfile, 
            file = sub(pattern = "\\.json$", replacement = ".csv", x = x))
})

jsonFile <- fromJSON(file="EA000408.json")
json_data_frame <- as.data.frame(jsonFile)
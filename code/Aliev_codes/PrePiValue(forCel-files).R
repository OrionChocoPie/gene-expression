# Проверка основных пакетов
if (!require("Biobase")) {
  source("http://bioconductor.org/biocLite.R");
  biocLite();
}
if (!require("affy")) {
  biocLite("affy");
}
if (!require("preprocessCore")) {
  biocLite("preprocessCore");
}
if (!require("gcrma")) {
  biocLite("gcrma");
}
if (!require("multtest")) {
  biocLite("multtest");
}
if (!require("stats")) {
  biocLite("stats");
}
if (!require("AnnotationDbi")) {
  biocLite("AnnotationDbi");
}

library("gcrma")
library("affyPrepr")


workingDir <- "C:/Users/Ruslan/Desktop/DiabetSets/11907/GSE11907_RAW"            #положение cel-файлов
if (!exists("workingDir"))
  stop("Define workingDir before sourcing this script.")

#names <- list.files(workingDir, pattern="*.CEL")                                 #alt описание sam, sbm в C:\Users\Ruslan\Desktop\DiabetSets\11907\...txt
#некоторые файлы имеют расширение .cel, а не .CEL

#names <- list.files(workingDir)

sam <-c(as.character(GSE11907.GPL96_series_matrix[219:370,1]),as.character(GSE11907.GPL96_series_matrix[218:369,2])) #sbm- вектор с названиями .CEL файлов для платформы (без '.CEL')
sam<-sort(sam)
sbm <-c(as.character(GSE11907.GPL97_series_matrix[188:308,1]),as.character(GSE11907.GPL97_series_matrix[187:307,2])) #sbm- вектор с названиями .CEL файлов для платформы (без '.CEL')
sbm<-sort(sbm)
#sbm<-paste0(sbm,'.CEL')
#names <- sbm
# Оказалось, что в файле описания больше названий, чем есть файлов .CEL на самом деле
names <- list.files(workingDir)
names <- substr(names,1,9)     #название файла .CEL до точки
names<-match(names,sbm)        #каждый элемент указывает на номер совпадения в sbm
names<-match(names,NA)
names<-match(names,NA)
names<-which (names==1)        #каждый элемент указывает на элементы в списке файлов рабочей папки, которые также есть и в sbm
names <- list.files(workingDir)[names]



files <- sapply(names, function (name) {
  file.path(workingDir, name)
})
t = load.affy(files)
                                                        #положение файла соответсвия
mapFile <- paste("C:/Users/Ruslan/Desktop/OncoFinder/packages/affyPrepr", "entrez_genesymbol_updated_2.txt", sep="/")
an = read.delim(mapFile, header=T, stringsAsFactors=F)
rownames(an) = an[, 1]
f = rpf(t, annotationTable=an, id="SYMBOL", unite="mean")

intensity = 2 ^ f

outFile <- file.path(workingDir, "11907B.txt")
write.table(intensity, outFile, sep="\t", quote=F)

print(paste("Created", outFile))

# Добавление атрибутов N и T из файла series_matrix.txt
Samples <- platform.HG.U133B
library("stringr")
normtum <- 'SYMBOL'
for (i in 1:(ncol(Samples)-1)){
  str1 <- GSE9006.GPL97_series_matrix[449+i%/%2,i%%2+1]
  print(str1)
  if (str_detect(str1,'Healthy'))
     normtum <- c(normtum,'Norm')
  else normtum <- c(normtum,'Tumour')}

Samples[1,] <- normtum
Samples <- Samples[,order('T' == substring(t(Samples[1,]),1,1))]
write.table(Samples,file=paste0('C:/Users/Ruslan/Desktop/DiabetSets/9006',"/9006B.txt"))

matr <- `9006A`
`9006A1` <- cbind(matr[,1],matr[,2:13],matr[,26:49])
`9006A2` <- cbind(matr[,1],matr[,14:25],matr[,50:72])
`9006A3` <- cbind(matr[,1],matr[,8:19],matr[,73:95])
`9006A4` <- cbind(matr[,1],matr[,2:7],matr[,20:25],matr[,96:118])
matr <- `9006B`
`9006B1` <- cbind(matr[,1],matr[,2:13],matr[,26:49])
`9006B2` <- cbind(matr[,1],matr[,14:25],matr[,50:72])
`9006B3` <- cbind(matr[,1],matr[,8:19],matr[,73:95])
`9006B4` <- cbind(matr[,1],matr[,2:7],matr[,20:25],matr[,96:118])


#Input
path <- 'C:/Users/Coolguy/Desktop/OncoFinder/datasets/1/'
name <- 'NB'
InSet <- read.table(paste0(path,name,'.txt'), quote="\"", stringsAsFactors=FALSE)

#QN all
a <- InSet
str <- a[1:nrow(a),1]
z<-ncol(a)
k <- which(substring(t(a[1,2:z]),1,1) == substring(t(a[1,2:z]),1,1)[z-1])[1] #Norm/Tumor column number
a <- apply(a[-1,-1],2,as.numeric)
a <- log10(a)
library(preprocessCore)
a <- normalize.quantiles(a[1:nrow(a),1:ncol(a)])
a <- 10^a 
a <- rbind('Tumour',a)
a <- cbind(as.character(str),a)
a[1,1:k]<-'Norm'
a[1,1]<-'SYMBOL'
write.table(a,file=paste0(path,name,"_QNall.txt"))

#QN Norm
a <- InSet
str <- a[1:nrow(a),1]
a <- apply(a[-1,-1],2,as.numeric)
a <- log10(a)
w <- 1
w[1:k] <- 1
w[(k+1):z] <- 0
library(preprocessCore)
a <- normalize.quantiles.robust(a[1:nrow(a),1:ncol(a)],w)
a <- 10^a 
a <- rbind('Tumour',a)
a <- cbind(as.character(str),a)
a[1,1:k]<-'Norm'
a[1,1]<-'SYMBOL'
write.table(a,file=paste0(path,name,"_QNnorm.txt"))

#QN Disease
a <- InSet
str <- a[1:nrow(a),1]
a <- apply(a[-1,-1],2,as.numeric)
a <- log10(a)
w <- 1
w[1:k] <- 0
w[(k+1):z] <- 1
library(preprocessCore)
a <- normalize.quantiles.robust(a[1:nrow(a),1:ncol(a)],w)
a <- 10^a
a <- rbind('Tumour',a)
a <- cbind(as.character(str),a)
a[1,1:k]<-'Norm'
a[1,1]<-'SYMBOL'
write.table(a,file=paste0(path,name,"_QNdisease.txt"))
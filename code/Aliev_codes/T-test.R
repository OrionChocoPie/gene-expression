#folder with datasets
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
paths <- 'C:/Users/Coolguy/Desktop/OncoFinder/datasets/1'
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


myfunc1 <- function(x,k,z){
  result <- (t.test(x[1:(k-1)],x[k:z])$p.value) #pi_value
  return(result)
}

myfunc2 <- function(x){
  result <- -log10(x)                           #-lg(pi_value)
  return(result)
}

myfunc3 <- function(x){
  result <- log10(x[1]/x[2])                    #lg(cnr)
  return(result)
}

myfunc4 <- function(x){
  result <- sign(x[3]-x[2])*x[1]                #-sign(T-N)*lg(pi_value)  x[1] from myfunc2
  return(result)
}

myfunc5 <- function(x,k,z){
  result <- (t.test(x[1:(k-1)],x[k:z])$p.value) #pi_value with CUT
  if (result > 0.05) result <- 1
  return(result)
}

path <- 'C:/Users/Anna_/Documents/GitHub/gene-expression/python/50'

#InSets <- list.files(paths)
#for (indsum in 1:length(InSets)){
#path <- paste0(paths,'/',substr(InSets[indsum],1,nchar(InSets[indsum])-4))

dir.create(path)
dir.create(paste0(path,'/log'))

InSet <- read.table(paste0(path,'.txt'), quote="\"", stringsAsFactors=FALSE,row.names=NULL)

InSet <- InSet[,-1]

InSet <- InSet[,order('T' == substring(t(InSet[1,]),1,1))]
# print(paste0(InSets[indsum],'    ',as.character(indsum),':',as.character(length(InSets))))
#empty space clean
InSet<-InSet[apply(InSet,1,sort)[1,]!='',]
#!!!!!!!!!!!Дебаг проверки на константную строку
InSet<-rbind(InSet,0)



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
  #pi_value for set`s pairs
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  maxnPairs <-400         #number of partitions for patients
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#time measure  
ptm <- proc.time()

#random sample generation
k <- which(substring(t(InSet[1,2:ncol(InSet)]),1,1) == substring(t(InSet[1,2:ncol(InSet)]),1,1)[ncol(InSet)-1])[1] #Norm/Tumor column number
nTum <- ncol(InSet)-k #number of 'Norm' 

#test for table`s size
if (choose(nTum,nTum%/%2)>(maxnPairs*2)) {nPairs <- maxnPairs} else {nPairs <- (choose(nTum,nTum%/%2)%/%2)+1}
n <- 1
vn <- ''
vn <- cbind(sort(sample(1:nTum, nTum%/%2, replace=F)))  #vector/number - vector contains random samples rooms
repeat{
  vn <- cbind(vn,sort(sample(1:nTum, nTum%/%2, replace=F)))
  fl <- FALSE
  for (i in 1:(ncol(vn)-1)) {if (all(vn[,i]==vn[,ncol(vn)])) fl <- TRUE}
  if (fl==TRUE) {vn <- vn[,1:(ncol(vn)-1)]
  print('wow')}
  else n <- n+1
  if (n%%25==0) print(n)
  if (n==nPairs) break
}
res<-''

#processing table pairs
for (n in 1:nPairs){
  vn1 <- c(1:k,k+vn[,n])
  a1 <- InSet[,vn1] #table1
  vn1 <- 1:nTum
  vn1 <- vn1[-vn[,n]]
  vn1 <- c(1:k,k+vn1)
  a2 <- InSet[,vn1] #table2
  
  #search/delete the constant vectors
  for (l in 1:2)
  {if (l == 1) a<-a1 else a<-a2 
  a <- a[c(TRUE,sapply(apply(a[-1,-1],1,unique),length) > 2),] #удаление константных строк
  z<-ncol(a)
  
  #processing
  a[1,z+1] <- 'pi_value'
  x<-apply(a[-1,-1],2,as.numeric)
  b <- apply(x,1,myfunc1,k,z)
  a[2:nrow(a),z+1] <- b
  
  a[1,z+2] <- '-lg(pi_value)'
  b <- apply(t(b),1,myfunc2)
  a[2:nrow(a),z+2] <- b
  
  a[1,z+3] <- 'mean(Norm)'
  x <- apply(a[-1,-c(1,(k+1):(z+3))],2,as.numeric)   #table with 'Norm' data 
  b <- apply(x,1,mean)
  a[2:nrow(a),z+3] <- b
  
  a[1,z+4] <- 'mean(Tumour)'  
  x <- apply(a[-1,-c(1:k,(z+1):(z+4))],2,as.numeric) #table with 'Tumour' data
  b <- apply(x,1,mean)
  a[2:nrow(a),z+4] <- b
  
  a[1,z+5] <- 'lg(cnr)'
  x <- apply(a[-1,(z+3):(z+4)],2,as.numeric)
  b <- apply(x,1,myfunc3)
  a[2:nrow(a),z+5] <- b
  
  a[1,z+6] <- '-sign(T-N)*lg(pi_value)'
  x <- apply(a[-1,(z+2):(z+4)],2,as.numeric)
  b <- apply(x,1,myfunc4)
  a[2:nrow(a),z+6] <- b
  
  if (l == 1) a1<-a else a2<-a}
  
  a <- merge(a1[2:nrow(a1),],a2[2:nrow(a2),], by = 1)  #!!!
  x <- as.numeric(a[,ncol(a1)])
  y <- as.numeric(a[,ncol(a)])
  res[n] <- cor(x, y, method = "pearson")
  if (n%%100==0) print(n)
}
write.table(vn,file=paste0(path,"/log/pv_vn.txt"))
write.table(res,file=paste0(path,"/log/pv_Result.txt"))
#создание densityplot
xlim<-as.numeric(c(min(as.numeric(res)),max(as.numeric(res))))
png(file=paste0(path,'/pi_value.png'))
res1 <- density(as.numeric(res),'bw' = 0.005) 
plot(res1,col='red',xlim=xlim)
dev.off()

print('pi_value')
print(proc.time() - ptm)



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#pvclust
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

#time measure
ptm <- proc.time()

library(pvclust)
library(magrittr)

k <- which(substring(t(InSet[1,2:ncol(InSet)]),1,1) == substring(t(InSet[1,2:ncol(InSet)]),1,1)[ncol(InSet)-1])[1] #Norm/Tumor column number
a<-InSet
colnames(a)<-substr(a[1,],3,length(a[1,]))
a <- a[-1,-(1:k)]
a <- apply(a,2,as.numeric)
a <- a[sapply(apply(a,1,unique),length) > 2,] #удаление константных строк
PVC <- pvclust(a, 'complete', 'euclidean', parallel = FALSE, nboot = 1000)
png(file=paste0(path,'/pvclust.png'))
#plot(PVC)
plot (PVC, hang = -1,main = "", sub = "", ylab = "Distance", xlab = "", col.pv=c("magenta","orange","grey"))
PVC %>% pvrect(border = "magenta")
dev.off()

print('pvclust')  
print(proc.time() - ptm)


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#just pi_value 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


a <- InSet
a <- a[c(TRUE,sapply(apply(a[-1,-1],1,unique),length) > 2),]
z<-ncol(a)
k <- which(substring(t(a[1,2:z]),1,1) == substring(t(a[1,2:z]),1,1)[z-1])[1] #Norm/Tumor column number

#processing
a[1,z+1] <- 'pi_value'
x<-apply(a[-1,-1],2,as.numeric)
b <- apply(x,1,myfunc1,k,z)
a[2:nrow(a),z+1] <- b

a[1,z+2] <- '-lg(pi_value)'
b <- apply(t(b),1,myfunc2)
a[2:nrow(a),z+2] <- b

a[1,z+3] <- 'mean(Norm)'
x <- apply(a[-1,-c(1,(k+1):(z+3))],2,as.numeric)   #table with 'Norm' data 
b <- apply(x,1,mean)
a[2:nrow(a),z+3] <- b

a[1,z+4] <- 'mean(Tumour)'  
x <- apply(a[-1,-c(1:k,(z+1):(z+4))],2,as.numeric) #table with 'Tumour' data
b <- apply(x,1,mean)
a[2:nrow(a),z+4] <- b

a[1,z+5] <- 'lg(cnr)'
x <- apply(a[-1,(z+3):(z+4)],2,as.numeric)
b <- apply(x,1,myfunc3)
a[2:nrow(a),z+5] <- b

a[1,z+6] <- '-sign(T-N)*lg(pi_value)'
x <- apply(a[-1,(z+2):(z+4)],2,as.numeric)
b <- apply(x,1,myfunc4)
a[2:nrow(a),z+6] <- b

write.table(a[,-(2:z)],file=paste0(path,"/log/just_pi_value.txt"))
#создание densityplot
xlim<-as.numeric(c(min(as.numeric(a[-1,ncol(a)])),max(as.numeric(a[-1,ncol(a)]))))
png(file=paste0(path,'/just_pi_value.png'))
res1 <- density(as.numeric(a[-1,ncol(a)]),'bw' = 0.005) 
plot(res1,col='red',xlim=xlim)
dev.off()
#}
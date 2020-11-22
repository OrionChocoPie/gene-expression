#folder with datasets
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
paths <- 'C:/Users/Coolguy/Desktop/OncoFinder/datasets/1'
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#library(stringr)
library(pvclust) #для pv_clust
library(magrittr) #для pv_clust
library(foreach) #для параллельного исчисления
library(doParallel) #для параллельного исчисления

#cores=detectCores()
cl <- makeCluster(4) 
registerDoParallel(cl)

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


InSets <- list.files(paths)
dir.create(paste0(paths,'/ZZZ_log_ZZZ'))
for (indsum in 1:length(InSets)){
  path <- paste0(paths,'/',substr(InSets[indsum],1,nchar(InSets[indsum])-4))
  InSet <- read.table(paste0(path,'.txt'), quote="\"", stringsAsFactors=FALSE,row.names=NULL)#row.names=NULL - защита от строк с одинаковым названием
  InSet <- InSet[,-1]#удаление столбца со старыми названиями строк
  InSet <- InSet[,order('T' == substring(t(InSet[1,]),1,1))]
  print(paste0(InSets[indsum],'    ',as.character(indsum),':',as.character(length(InSets))))
  #empty space clean
  InSet<-InSet[apply(InSet,1,sort)[1,]!='',]
  #!!!!!!!!!!!Дебаг проверки на константную строку
  InSet<-rbind(InSet,0)
  
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #pi_value for set`s pairs
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  maxnPairs <-200         #number of partitions for patients
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #random sample generation
  k <- which(substring(t(InSet[1,2:ncol(InSet)]),1,1) == substring(t(InSet[1,2:ncol(InSet)]),1,1)[ncol(InSet)-1])[1] #Norm/Tumor column number
  nTum <- ncol(InSet)-k #number of 'Norm' 
  
  #test for table`s size
  if (choose(nTum,nTum%/%3)>(maxnPairs*2)) {nPairs <- maxnPairs} else {nPairs <- (choose(nTum,nTum%/%3)%/%2)+1}
  n <- 1
  vn <- ''
  vn <- cbind(sort(sample(1:nTum, nTum%/%3, replace=F)))  #vector/number - vector contains random samples rooms
  repeat{
    vn <- cbind(vn,sort(sample(1:nTum, nTum%/%3, replace=F)))
    fl <- FALSE
    for (i in 1:(ncol(vn)-1)) {if (all(vn[,i]==vn[,ncol(vn)])) fl <- TRUE}
    if (fl==TRUE) vn <- vn[,1:(ncol(vn)-1)]
    else n <- n+1
    if (n==nPairs) break
  }
  res<-''
  
  #processing table pairs
  res <- foreach (n = 1:nPairs, .combine=cbind) %dopar% {
    vn1 <- c(1:k,k+vn[,n])
    a1 <- cbind(InSet[,vn1],InSet[,k+vn[,n]]) #table1, doubled
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
    cor(x, y, method = "pearson")
  }
  #создание densityplot
  xlim<-as.numeric(c(min(as.numeric(res)),max(as.numeric(res))))
  png(file=paste0(paths,'/ZZZ_log_ZZZ/',substr(InSets[indsum],1,nchar(InSets[indsum])-4),'.png'))
  res1 <- density(as.numeric(res),'bw' = 0.005) 
  plot(res1,col='red',xlim=xlim)
  dev.off()
  

  
}
stopCluster(cl)

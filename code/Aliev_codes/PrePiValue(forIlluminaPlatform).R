#С таблицей перевода в гены
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
a<-read.table('C:/Users/Coolguy/Desktop/OncoFinder/datasets/Venous_thromboembolism/V_Thr_Set.txt', quote="\"", stringsAsFactors=FALSE)
a <- GSE30210_series_matrix_short
SYMBOLS <- SYMBOLS
Sample <- SampleInfo
Path <- 'C:/Users/Ruslan/Desktop/OncoFinder/datasets/DiabetSets/30210/30210.txt'
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

library("stringr")

#добавление названия генов

#x<-merge(SYMBOLS,cbind(a[2:nrow(a),1]), by ='V1')
a<-merge(SYMBOLS,a, by ='V1')

#добавление Norm/Tumour
for (i in 1:(ncol(a)-2)){
  if (str_detect(Sample[i,2],'Control'))
    a[1,2+i] <- 'Norm'
  else a[1,2+i] <- 'Tumour'}
a <- a[,order('T' == substring(t(a[1,]),1,1))]
#удаление столбца ID_REF
a <- a[,2:ncol(a)]
colnames(a)<-'V1'
write.table(a,file=Path)



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Без таблицы перевода в гены (N/T)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
a<-read.table('C:/Users/Coolguy/Desktop/OncoFinder/datasets/Microglia/GSE1432_series_matrix.txt', quote="\"", stringsAsFactors=FALSE)
Sample<-read.table('C:/Users/Coolguy/Desktop/OncoFinder/datasets/Microglia/SetInfo.txt', quote="\"", stringsAsFactors=FALSE)
#
#a <- `Thr`
#Sample <- `SetInfo`
#
Path <- 'C:/Users/Coolguy/Desktop/OncoFinder/datasets/Microglia/Microglia.txt'
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library("stringr")
#добавление Norm/Tumour
for (i in 1:(ncol(a)-1)){
  if (str_detect(Sample[i,2],'Control')) #Control, control, non disease name, healthy, space+control, Ctr 
    a[1,1+i] <- 'N_GSE1432'
  else a[1,1+i] <- 'T_GSE1432'}
a[1,1]<-'SYMBOL'
a <- a[,order('T' == substring(t(a[1,]),1,1))]
write.table(a,file=Path)


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Без таблицы перевода в гены (больше состояний)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Path <- 'C:/Users/Coolguy/Desktop/OncoFinder/datasets/Venous_thromboembolism/GSE48000/'
a <- read.table(paste0(Path,'/GSE48000.txt'), quote="\"", stringsAsFactors=FALSE)
Sample <- read.table(paste0(Path,'/SetInfo.txt'), quote="\"", stringsAsFactors=FALSE)
name <- 'Thr_GSE48000'
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library("stringr")
#добавление Norm/Tumour
a[1,1]<-'SYMBOL'
`High` <- a[,1]
`Low` <- a[,1]
`Moderate` <- a[,1]
for (i in 1:(ncol(a)-1)){
  if (str_detect(Sample[i,2],'Healthy')){
    a[1,1+i] <- 'Norm'
    `High` <- cbind(High,a[,1+i])
    `Low` <- cbind(`Low`,a[,1+i])
    `Moderate` <- cbind(Moderate,a[,1+i])}
  if (str_detect(Sample[i,2],'High')){
    a[1,1+i] <- 'Tumour_High'
    High <- cbind(High,a[,1+i])}
  if (str_detect(Sample[i,2],'Low')){
    a[1,1+i] <- 'Tumour_Low'
    `Low` <- cbind(`Low`,a[,1+i])}
  if (str_detect(Sample[i,2],'Moderate')){
    a[1,1+i] <- 'Tumour_Moderate'
    Moderate <- cbind(Moderate,a[,1+i])}}
#сортировка и сохранение
a <- a[,order('T' == substring(t(a[1,]),1,1))]
High <- High[,order('T' == substring(t(High[1,]),1,1))]
`Low` <- `Low`[,order('T' == substring(t(`Low`[1,]),1,1))]
Moderate <- Moderate[,order('T' == substring(t(Moderate[1,]),1,1))]
colnames(High)[1] <- colnames(`Low`)[1] <- colnames(Moderate)[1] <- 'V1'
write.table(a,file=paste0(Path,name,'_all.txt'))
write.table(High,file=paste0(Path,name,'_High.txt'))
write.table(`Low`,file=paste0(Path,name,'_Low.txt'))
write.table(Moderate,file=paste0(Path,name,'_Moderate.txt'))

#<0 clear
path <- 'C:/Users/Coolguy/Desktop/OncoFinder/datasets/1/'
name <- 'NB'
#InSet <- read.table(paste0(path,name,'.txt'), quote="\"", stringsAsFactors=FALSE)
InSet <- NB
a <- InSet[-1,]
a <- a[apply(apply(a[,-1],2,as.numeric),1,sort)[1,]>0,]
a <- rbind(InSet[1,],a)
write.table(a,file=paste0(path,name,".txt"))


#На заметку, как убрать цикл для определения Norm/Tumour (без отдельного файла, если с ним, то надо t(str_detect(...)))
library("stringr")
a[1,str_detect(a[1,],'nor')]<-paste0('N_',a[1,str_detect(a[1,],'nor')])  #Control, control, non disease name, healthy, space+control, Ctr, nor
a[1,!(str_detect(a[1,],'nor'))]<-paste0('T_',a[1,!(str_detect(a[1,],'nor'))])

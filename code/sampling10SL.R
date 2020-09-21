library(e1071)
library(rpart)
library(caTools)
library(preprocessCore)
library(stats)
#library(rJava)
#library(xlsxjars)
#library(xlsx)
library(rjson)
library(XML)
library(plyr)
library(stringr)
library(dendextend)
library(pvclust)
library(oligo)
library(stringr)
library(hgu133plus2hsentrezgprobe)
library(pd.hgu133plus2.hs.entrezg)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
library(affy)
library(hugene11stv1cdf)

#prefix = '/media/nicolas/795273AB6F22B1C1/New_archive/cancer_harmony'
prefix = "/mnt/disk1/New_archive/Sh_SEQC_GDC"
#prefix = "/home/nborisov/gdc"

Ps = c('571', '10558', '11154', 'AGL', '75')
NP = length(Ps)

Qs = c("0", "1")
NQ = length(Qs)

for ( nq in 1:NQ ) {

   for ( np in 1:NP ) {
   
        IFN = paste(prefix, "/Sh_", Ps[np], "_", Qs[nq], ".csv", sep = "")
        Sh = read.table(IFN, header = TRUE, sep = ",")
        
        print (np)

        for ( iii in 1:25 ) {
        
             print (iii)
                 
             if ( iii < 10 )  IFN = paste(prefix, "/reducedSL/meta_0", toString(iii), ".csv", sep = "")
             if ( iii >= 10 )  IFN = paste(prefix, "/reducedSL/meta_", toString(iii), ".csv", sep = "")
             
             meta = read.table(IFN, header = TRUE, sep = ",")
             
             SN = as.vector(meta[,1])
             
             NS = length(SN)
             
             for ( ns in 1:NS ) {
             
                 L = str_length(SN[ns])
                 STR = str_sub(SN[ns], start = 3, end = L)
                 NUM = as.numeric(STR)
                 
                 if ( ns == 1 ) INDEX = NUM
                 if ( ns > 1 ) INDEX = c(INDEX,NUM)

             }
             
             sh = Sh[,INDEX]


             if ( iii < 10 )  { 
     
                  OFN = paste(prefix, "/reducedSL/Sh_", Ps[np], "_", Qs[nq], "_0", toString(iii), ".csv", sep = "")
                  write.table(sh, OFN, col.names = SN, row.names = TRUE, sep = ",")
             }
    
             if ( iii >= 10 )  { 
                  OFN = paste(prefix, "/reducedSL/Sh_", Ps[np], "_", Qs[nq], "_", toString(iii), ".csv", sep = "")
                  write.table(sh, OFN, col.names = SN, row.names = TRUE, sep = ",")        
             }
             
       }      
       
   }    

}



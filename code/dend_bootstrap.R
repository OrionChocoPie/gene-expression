#library(e1071)
library(rpart)
#library(caTools)
#library(preprocessCore)
library(stats)
#library(rJava)
#library(xlsxjars)
library(openxlsx)
library(ggplot2)
library(vioplot)
library(dendextend)
library(stringr)
library(pvclust)
library(stats)
#library(DESeq)
#library(DESeq2)
#library(oligo)
#library(CONOR)
#library(HARMONY)

bar_y_scale = function(dendrogram) {
    return(attr(dendrogram, "height") * 0.1)
}

prefix = "/mnt/disk1/New_archive/paired"

Ps = c('0', '1' '2', '3', '4', '5', '6')
NP = length(Ps)

Qs = c("0", "1")
NQ = length(Qs)

methods = array("XYZXYZ", dim = c(NP*NQ))


for ( nq in 1:NQ ) {

   for ( np in 1:NP ) {
   
        methods[np + (nq-1)*NP] = paste(Ps[np], Qs[nq], sep = "_")
        
   }

}

M = length(methods)

bar_row_labels = c("Platform", "Samples")

title_for_method = c("Sh20","Sh30","Sh40","Sh50","Sh60","Sh21","Sh31","Sh41","Sh51","Sh61") 

#for ( iii in 1:25 ) {
iii = 1

if (iii < 10 ) IFN = paste(prefix, "/reduced/meta_0", toString(iii), ".csv", sep = "")
if (iii >= 10 ) IFN = paste(prefix, "/reduced/meta_", toString(iii), ".csv", sep = "")
meta = read.table(IFN, header = TRUE, sep = ",")

NS = nrow(meta)

D = as.vector(meta[,2])
P = as.vector(meta[,3])
S = as.vector(meta[,1])

UD1 = unique(D)
UP1 = unique(P)

keys = c(8,10,11,9,1,2,3,4,5,6,7)
UD = UD1[keys]

UP = sort(UP1)

NUD = length(UD)
CT_D = 1:NUD
t_rgb = col2rgb("red")
t_hsv = rgb2hsv(t_rgb)     
for ( nud in 1:NUD ) {
     if (nud <= 4) t_hsv[1] = (nud-1)/(NUD+4);
     if (nud >= 5) t_hsv[1] = (nud+2)/(NUD+4);
     CT_D[nud] = hsv(t_hsv[1],1,t_hsv[3]);
}

NUP = length(UP)
CT_P = 1:NUP
t_rgb = col2rgb("red")
t_hsv = rgb2hsv(t_rgb)     
for ( nup in 1:NUP ) {
     t_hsv[1] = (nup-1)/(NUP+1);
     CT_P[nup] = hsv(t_hsv[1],0.33,t_hsv[3]);
}

CD = array(dim=c(NS))
CP = array(dim=c(NS))

for ( ns in 1:NS ) {
     d = D[ns]
     index = which (d == UD)
     CD[ns] = CT_D[index]
     p = P[ns]
     index = which (p == UP)
     CP[ns] = CT_P[index]
}

for ( m in 1:M )  {

  print (m)

  if (iii < 10 ) IFN = paste(prefix, "/reduced/", methods[m], "_0", toString(iii), ".csv", sep = "")
  if (iii >= 10 ) IFN = paste(prefix, "/reduced/", methods[m], "_", toString(iii), ".csv", sep = "")
  mas_all = read.table(IFN, header = TRUE, sep = ",")
  
  mas_all = as.matrix(mas_all) 
  mas_all = log10(mas_all+1)
  mas_all = t(mas_all)

  DDD = dist(mas_all)
  par(cex.axis = 1) 
  hc = hclust(DDD)
  hcd = as.dendrogram (hc)
  dend_labels = hcd %>% labels
  nlabels = length(dend_labels)
  bar_color_cancer = array("XYZXYZ", dim = c(nlabels))
  bar_color_platform = array("XYZXYZ", dim = c(nlabels))
  for ( ilabel in 1:nlabels ) {
       sample = dend_labels[ilabel]
       index = which ( sample == S )
       bar_color_cancer[ilabel] = CD[index]
       bar_color_platform[ilabel] = CP[index]
       dend_labels[ilabel]=""
  }     
  if ( iii < 10 ) GFN = paste(prefix, "/dend/", methods[m], "_0", toString(iii),  ".ps", sep = "")
  if ( iii >= 10 ) GFN = paste(prefix, "/dend/", methods[m], "_", toString(iii),  ".ps", sep = "")
  postscript(file=GFN)
  par(cex.axis = 1.5) 
  par(cex.main = 2)
  hcd %>% set("labels", dend_labels) %>% plot()
  par(cex.lab = 1.5)
  par(mar = c(8,4,4,5)+.1)
  title(main = title_for_method[m], sub = "", ylab = "Distance")

  bar_colors = cbind(bar_color_platform, bar_color_cancer)
  colored_bars(
        colors=bar_colors,
        dend=hcd,
        sort_by_labels_order = FALSE,
        rowLabels = bar_row_labels,
        cex.rowLabels = 1.5,
        y_scale = bar_y_scale(hcd),
        y_shift = -bar_y_scale(hcd)
  )
  dev.off()

  PVC <- pvclust(t(mas_all), 'complete', 'euclidean', parallel = TRUE, nboot = 300, quiet=FALSE)
  
  hcd = as.dendrogram (PVC)
  dend_labels = hcd %>% labels
  
   nlabels = length(dend_labels)
   dend_labels = hcd %>% labels
   nlabels = length(dend_labels)
   bar_color_cancer = array("XYZXYZ", dim = c(nlabels))
   bar_color_platform = array("XYZXYZ", dim = c(nlabels))
   for ( ilabel in 1:nlabels ) {
       sample = dend_labels[ilabel]
       index = which ( sample == S )
       bar_color_cancer[ilabel] = CD[index]
       bar_color_platform[ilabel] = CP[index]
       dend_labels[ilabel]=""
   } 
   
   if ( iii < 10 ) GFN = paste(prefix, "/dend/", methods[m], "_0", toString(iii),  "_rect.ps", sep = "")
   if ( iii >= 10 ) GFN = paste(prefix, "/dend/", methods[m], "_", toString(iii),  "_rect.ps", sep = "")
   postscript(file=GFN)
   postscript(file=GFN)
   par(cex.axis = 1.5) 
   par(cex.main = 2)
   hcd %>% set("labels", dend_labels) %>% plot()
   par(cex.lab = 1.5)
   par(mar = c(8,4,4,5)+.1)
   title(main = title_for_method[m], sub = "", ylab = "Distance")

   text (PVC,col = c("darkolivegreen4","orange","grey"))
   PVC %>% pvrect(border = "magenta")
   
  bar_colors = cbind(bar_color_platform, bar_color_cancer)
  colored_bars(
        colors=bar_colors,
        dend=hcd,
        sort_by_labels_order = FALSE,
        rowLabels = bar_row_labels,
        cex.rowLabels = 1.5,
        y_scale = bar_y_scale(hcd),
        y_shift = -bar_y_scale(hcd)
  )
  dev.off()

}  

#}

GFN = paste(prefix, "/dend/legend10M.ps", sep = "")
postscript(file=GFN)
plot.new( )
legend("topleft", legend = UD, fill = CT_D, border = CT_D, bty = "n", cex = 1.9)
legend("topright", legend = UP, fill = CT_P, border = CT_P, bty = "n", cex = 2) 
dev.off()
  

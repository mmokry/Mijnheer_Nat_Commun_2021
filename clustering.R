library(rpart)
library(scatterplot3d)
library(limma)
library(gplots)
library(IDPmisc)


list = read.table (file = "list.txt",row.names=1) #list with individual RPKM files, one tab delimited file per sample.  each row is different genes, columns are as follow: "gene","RS ID","NA","chr","start","end","strand","transcript_length","RPKM"


a= read.table (file = row.names(list)[1],sep = "\t",header = F)
len = length(row.names(list))

for (i in 1:len){
  print(row.names(list)[i])
  a[,(i+8)]=read.table (file = row.names(list)[i],sep = "\t",header = F)[9]
}


list = read.table (file = "names.txt",row.names = 1,sep = "\t")
names(a)= c("gene","RS","NA","chr","start","end","strand","transcript_length",row.names(list))


colors= read.table (file = "col.txt") #color encoding for graphs per sample. Same order as "list.txt", each row is one sample, colors are encoded as "1", "2". etc

#quantile normalizations
ln = length(a[1,])
a[,9:ln]=normalizeQuantiles(a[,9:ln])

#log2 scaling
a[,9:ln] = log2(a[,9:ln]+0.1)   #log2 scale for RPKM
b=data.matrix(a[,9:ln])

#save normalize file for the reference
write.table(a,file = "GM_RPKM_QN_rna.txt", quote = F, row.names = F, sep =  "\t" )

#select variable genes
    fc = 1
    min_val = 2
    uq =  0.90   # upper quantile treshold
    lq = 0.10     #lower quantile treshold
    m = 0
    count = 0
    bbb=b
    #check if there is outlier
    for (i in 1:length(a[,1])){
      if (quantile(b[i,],c(uq))>(fc+(quantile(b[i,],c(lq))))){
        if ((quantile(b[i,],c(uq))>min_val)){
#	if (max(b[i,]>min_val)){
	 count = count+1   
	 m[count]=i
	}
      }
    b[i,]=b[i,]-median(b[i,])
    bbb[i,1]=median(b[i,1:3]) #it meas that samples 1:3 are in one group
   bbb[i,2]=median(b[i,4:7]) 
   bbb[i,3]=median(b[i,8:11])

    }
print(length(m))
################################### heatmaps

pdf(file = "heatmap2_cor.pdf",width =7, height = 7)
heatmap.2(cor(b[m,],method = "pearson"),trace = "none",margins = c(10,10),density.info = "none",col = colorpanel(1024, "blue","black","yellow"))
dev.off()

pdf(file = "heatmap2_genes.pdf",width =7, height = 10)
heatmap.2((b[m,]),trace = "none",margins = c(10,10),density.info = "none",scale = "row",col = colorpanel(1024, "blue","black","yellow"))
dev.off()
####################################PCA

pdf(file = "PCA_3D.pdf", width = 5, height = 5)
PCA = prcomp(t(b[m,1:(ln-8)]))

PC1L=paste("PC1 (",round(summary(PCA)$importance[2,1]*100,digits = 1),"%)",sep = ""  )
PC2L=paste("PC2 (",round(summary(PCA)$importance[2,2]*100,digits = 1),"%)",sep = ""  )
PC3L=paste("PC3 (",round(summary(PCA)$importance[2,3]*100,digits = 1),"%)",sep = ""  )


RR = scatterplot3d(PCA$x[,1:3], type = "h",xlab = PC1L,ylab = PC2L,zlab = PC3L)


for (i in 1:(ln-8)){
  x = (data.matrix(PCA$x)[,1])
  y = (data.matrix(PCA$x)[,2])
  z = (data.matrix(PCA$x)[,3])
  RR$points3d(RR$points3d(data.matrix(x[i]),data.matrix(y[i]),data.matrix(z[i]),
	    pch = 20, col = colors[i,1]), pch = 20, col = colors[i,1])
  print(colors[i,1])
}
dev.off()


pdf(file = "PCA_1_2_unsupervised.pdf", width = 7, height = 7)
plot(PCA$x[,c(1,2)],col = colors[,1],pch = 20,cex = 2.5,xlab = PC1L,ylab = PC2L);legend("topleft",
legend=c("childPBTreg","JIASFTeff","JIASFTreg"),col = c(4:6),pch = 20)
points(PCA$x[,c(1,2)],col = 1,pch = 1,cex = 2)
dev.off()

pdf(file = "PCA_1_3_unsupervised.pdf", width = 7, height = 7)
plot(PCA$x[,c(1,3)],col = colors[,1],pch = 20,cex = 2.5,xlab = PC1L,ylab = PC3L);legend("bottomright",
legend=c("childPBTreg","JIASFTeff","JIASFTreg"),col = c(4:6),pch = 20)
points(PCA$x[,c(1,3)],col = 1,pch = 1,cex = 2)
dev.off()

pdf(file = "PCA_2_3_unsupervised.pdf", width = 7, height = 7)
plot(PCA$x[,c(2,3)],col = colors[,1],pch = 20,cex = 2.5,xlab = PC2L,ylab = PC3L);legend("bottomright",
legend=c("childPBTreg","JIASFTeff","JIASFTreg"),col = c(4:6),pch = 20)
points(PCA$x[,c(2,3)],col = 1,pch = 1,cex = 2)
dev.off()

############################ k-means

cl = 14; #number of groups

b=b[m,]
bbb=bbb[m,]
set.seed(10)
r =  kmeans(bbb[,1:3],centers = cl,nstart = 2000)

pdf(file = "Heatmap_Clusters_all.pdf",width = 4, height = 8)
par(mar = c(0,0,0,0))
image(t(b[order(r$cluster),]),col = colorpanel(100, "blue", "black","yellow"),zlim = c(min(b[,])/1.5,max(b[,])/1.5))
dev.off()
a = a[m,]
a[,9:ln]=b
a[,(ln+1)] =r$cluster

b=data.matrix(a[,9:(ln+1)])

ln = ln -8


for (i in 1:cl){
 
 pdf(file = paste("Boxplot_Cluster_",i,".pdf",sep = ""),width = 5, height = 5)
  boxplot(b[,1][b[,(ln+1)]==i],b[,2][b[,(ln+1)]==i],b[,3][b[,(ln+1)]==i],b[,4][b[,(ln+1)]==i],  b[,5][b[,(ln+1)]==i],b[,6][b[,(ln+1)]==i],b[,7][b[,(ln+1)]==i],b[,8][b[,(ln+1)]==i], b[,9][b[,(ln+1)]==i],b[,10][b[,(ln+1)]==i],b[,11][b[,(ln+1)]==i] ,  
   ylim = c(min(b[,1:ln]),max(b[,1:ln])),ylab = "log2RPKM (median centered)",main = paste("Cluster",i,"\n","n=",length(b[,9][b[,(ln+1)]==i]),sep = ""),col = c(4,4,4,5,5,5,5,6,6,6,6))
  dev.off()


  pdf(file = paste("Heatmap_Cluster_",i,".pdf",sep = ""),width = 5, height = length(b[,12][b[,(ln+1)]==i])/1000)
  par(mar = c(0,0,0,0))
  image(t(array(b[,1:ln][b[,(ln+1)]==i],dim = c(length(b[,1:ln][b[,(ln+1)]==i])/ln,ln))),col = colorpanel(100, "blue", "black","yellow"),zlim = c(min(b[,1:ln])/1.5,max(b[,1:ln])/1.5))
  dev.off()
  print (paste(i," ",length(b[,9][b[,(ln+1)]==i])))
  
  
}
pdf(file = "Boxplot_Clusters_all.pdf",width = 15, height = 12)
par (mfrow = c(4,5))
for (i in 1:cl){
   boxplot(b[,1][b[,(ln+1)]==i],b[,2][b[,(ln+1)]==i],b[,3][b[,(ln+1)]==i],b[,4][b[,(ln+1)]==i],  b[,5][b[,(ln+1)]==i],b[,6][b[,(ln+1)]==i],b[,7][b[,(ln+1)]==i],b[,8][b[,(ln+1)]==i], b[,9][b[,(ln+1)]==i],b[,10][b[,(ln+1)]==i],b[,11][b[,(ln+1)]==i] ,  
   ylim = c(min(b[,1:ln]),max(b[,1:ln])),ylab = "log2RPKM (median centered)",main = paste("Cluster",i,"\n","n=",length(b[,9][b[,(ln+1)]==i]),sep = ""),col = c(4,4,4,5,5,5,5,6,6,6,6))
  }

dev.off()
###########################writecluster file
list = read.table (file = "list.txt",row.names=1)

a= read.table (file = row.names(list)[1],sep = "\t",header = F)
len = length(row.names(list))

for (i in 1:len){
  print(row.names(list)[i])
  a[,(i+8)]=read.table (file = row.names(list)[i],sep = "\t",header = F)[9]
}


list = read.table (file = "names.txt",row.names = 1,sep = "\t")
names(a)= c("gene","RS","NA","chr","start","end","strand","transcript_length",row.names(list))
ln = length(a[1,])
a[,9:ln]=normalizeQuantiles(a[,9:ln])

a[,9:ln] = log2(a[,9:ln]+0.1)   #log2 scale for RPKM

a = a[m,]
a[,(ln+1)]=r$cluster

write.table(a,file ="Genes_Clusters.txt",sep = "\t",quote = F, col.names= T  )


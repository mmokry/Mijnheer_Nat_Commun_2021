library(DESeq2)
a = read.table (file = "Raw.counts_rna.txt",header = T,sep = "\t",row.names = 1) #tab delimited gene count table file with samples in columns and genes in rows

##########


coldata = read.table(file = "names_DS.txt",header = T,sep = "\t")[c(12:14,19:22),c(1,3,2)] #annotation file with samples as rows and parameters/groups as columns (only samples used in the analysis are selected)

b=a[,c(19:21,26:29)] #subselection of samples from gene count table

dso = DESeqDataSetFromMatrix(countData = b,colData = coldata,design =  ~pos)
dsoA<-DESeq(dso)

#MA_plot
pdf(file = "A_vs_B_MAplot.pdf")
res = results(dsoA)
res$pval[is.na(res$pval)]<-1
res$padj[is.na(res$padj)]<-1
plot(log2(results(dsoA)$baseMean),results(dsoA)$log2FoldChange,pch = 20,cex = 0.3,col=ifelse(res$pval>=0.05,"black","red"),ylim = c(-8,8),xlab = "log2 baseMean",ylab = "log2 fold change");abline(h=0, col = "red",lwd = 2)
points(log2(results(dsoA)$baseMean),results(dsoA)$log2FoldChange,pch = 20,cex = ifelse(res$baseMean<=10,0.4,0),col=ifelse(res$baseMean<=10,"black","red"))
points(log2(results(dsoA)$baseMean),results(dsoA)$log2FoldChange,pch = 20,cex = ifelse(abs(res$log2FoldChange)<=0.6,0.4,0),col=ifelse(abs(res$log2FoldChange)<=0.6,"black","red"))
dev.off()



#volcano_plot
pdf(file = "1_A_vs_B_Vplot.pdf",width =5, height = 5)
plot(results(dsoA)$log2FoldChange,log10(1/(res$pval)),pch = 20,cex = 0.4,col=ifelse(res$pval>=0.05,"black","red"),ylim = c(0,max(log10(1/(res$pval)))),xlab = "log2 fold change",ylab = "-log10(pval)");abline(v=0, col = "grey",lwd = 2);abline(h=log10(1/0.05),col = "red")
points(results(dsoA)$log2FoldChange,log10(1/(res$pval)),pch = 20,cex = ifelse(res$baseMean<=10,0.4,0),col=ifelse(res$baseMean<=10,"black","red"),ylim = c(0,max(log10(1/(res$pval)))))
points(results(dsoA)$log2FoldChange,log10(1/(res$pval)),pch = 20,cex = ifelse(abs(res$log2FoldChange)<=0.6,0.4,0),col=ifelse(abs(res$log2FoldChange)<=0.6,"black","red"),ylim = c(0,max(log10(1/(res$pval)))))
dev.off()

#PCA_plot
pdf(file = "3_childPBTreg_vs_JIASFTreg_PCAplot5000.pdf")
plotPCA(rlog(dsoA),intgroup=c("pos"),ntop = 5000)
dev.off()

#output_DEGs
write.table((results(dsoA)),file="1_A_vs_B.txt",sep = "\t",quote = F)

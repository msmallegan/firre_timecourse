# written by Taeyoung Hwang, taeyoungh@gmail.com
# last update: 5/1/2018

fCountReader <- function (countDir, sampleID, suffix=".fCounts") {
  
  for (i in 1:length(sampleID)) {
    filename=paste(countDir,"/",sampleID[i],suffix,sep="")
    print(filename)
    temp=read.table(filename,header=T,stringsAsFactors = F)
    if (i==1) {
      annot=temp[,c(1,5,6)]
      fCount=data.frame(temp[,8])
    } else {
      fCount=data.frame(fCount,temp[,8])
    }
  }
  colnames(annot)=c("geneID","strand","length")
  annot$strand=substr(annot$strand,1,1)
  
  colnames(fCount)=sampleID
  rownames(fCount)=annot$geneID
  return(list(count=fCount,annot=annot))
}

countConverter<-function(fCount,return="TPM") {
  if (return=="FPKM") {
    out=fCount$count/(fCount$annot$length/1000)
    out=t(t(out)/(colSums(fCount$count)/10^6))
  } else if (return=="TPM") {
    out=fCount$count/(fCount$annot$length/1000)
    out=t(t(out)/(colSums(out)/10^6))
  }
  return(out)
}

plotPCA <- function (expr, sampleSheet, colorVar=NA, shapeVar=NA, pointSize=5) {
  require("ggplot2")
  pcaObj<-prcomp(t(expr[which(apply(expr,1,var)>0),]), center=T, scale=T)
  pcaRes=data.frame(sampleSheet,PC1=pcaObj$x[,1],PC2=pcaObj$x[,2])
  pcaVarsPct=signif(((pcaObj$sdev)^2)/(sum((pcaObj$sdev)^2)),3)*100
  
  if (!is.na(colorVar) & !is.na(shapeVar)) {
    p=ggplot(pcaRes, aes_string("PC1", "PC2", color=colorVar, shape=shapeVar))
  }  
  
  if (!is.na(colorVar) & is.na(shapeVar)) {
    p=ggplot(pcaRes, aes_string("PC1", "PC2", color=colorVar))
  }  
  
  if (is.na(colorVar) & !is.na(shapeVar)) {
    p=ggplot(pcaRes, aes_string("PC1", "PC2", shape=shapeVar))
  }  
  
  if (is.na(colorVar) & is.na(shapeVar)) {
    p=ggplot(pcaRes, aes_string("PC1", "PC2"))
  }  
  
  p + geom_point(size=pointSize) +
    xlab(paste0("PC1: ",pcaVarsPct[1],"% variance")) + ylab(paste0("PC2: ",pcaVarsPct[2],"% variance"))
  
}


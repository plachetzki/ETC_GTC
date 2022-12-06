

#source("https://bioconductor.org/biocLite.R"); biocLite("tximport"); install.packages("readr"); biocLite("edgeR")    #<-----uncomment & run if any not installed yet
#Bioconductor our out of date- run code below to install newer version 
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install(version = "3.12")
#BiocManager::install(c("tximport", "edgeR", "biocLite"))
library(tidyr); library(dplyr)
library(cowplot)
library(ape)
library(ggtree)
library(tximport); library(readr); library(edgeR); library(limma)
library(ggplot2)
library(this.path)
dir <-this.dir()
setwd(dir)
list.files()

#in your wd make a text file named libraries_to_stages.txt - this should have 2 columns - the first is the name
#of your sample and the second (separated by a space) is the group number it belongs to (so the stage)
#Make a dir in your cwd named mapping - place all of your output dirs for each sample that was outputted by salmon in this dir

#library info file
sample_meta<-read.table("libraries_to_stages_EG_ind4.txt",header=F,row.names=1)
colnames(sample_meta)<-"group"
skin<-rownames(sample_meta)[which(sample_meta$group=="skin")]; skinfiles <- file.path(dir, "mapping",skin, "quant.sf"); names(skinfiles)<-skin
sg<-rownames(sample_meta)[which(sample_meta$group=="slime_gland")]; sgfiles <- file.path(dir, "mapping",sg, "quant.sf"); names(sgfiles)<-sg

pdffile<-"EG_ind4.pdf"
##pick colors for each library type
skincol <- "darkgreen"
sgcol <- "darkred"


quantfiles<-c(skinfiles,sgfiles); 
cols<-rep(c(skincol, sgcol),each=3) #change the last number to how many groups you have

################################
### Multidimensional scaling ###
################################

runEdgeRmds<-function(salmonquantfiles, colors,meta){
  ##  read in files with tximport
  txi.salmon<- tximport(salmonquantfiles, type = "salmon", txOut=T, dropInfReps=TRUE)
  cts <- txi.salmon$counts
  print(colSums(cts))
  
  keep <- rowSums(cpm(cts)>10.0) >= nrow(meta)
  cts<- cts[keep,]
  dim(cts)
  print(colSums(cts))
  
  group <- groups
  y <- DGEList(counts=cts ,group=meta$group)
  y <- calcNormFactors(y)
  y <-estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y, prior.n=16) # TagwiseDisp n-value should be close to: 50/(#samples - #groups) = 50/(36-6) = 50/30 =1.66
  
  
  plotMDS.DGEList(y , main = "MDS Plot for Count Data", labels = colnames(y),col=colors, xpd = TRUE)
  legend(1.25,0.2 ,legend=paste("Stage ",levels(as.factor(meta$group))),text.col=colors[seq(1,length(colors),3)], xpd = TRUE, cex = 0.75) #change the last number to how many groups you have
  
}

#runEdgeRmds(quantfiles,cols,sample_meta)
##############################################
### pairwise comparisons between libraries ###
##############################################

runEdgeRpwise<-function(salmonquantfiles,colors,meta){
  ##  read in files with tximport
  txi.salmon<- tximport(salmonquantfiles, type = "salmon", txOut=T, dropInfReps=TRUE)
  cts <- txi.salmon$counts
  print(colSums(cts))
 # tags2find<-read.csv("EGtags.csv",header=T) #MSP 
 # cts[tags2find$transcript.name,]
  #keep <- rowSums(cpm(cts)>10.0) >= length(groups) #keeps tags which are 10+ in every library
  keep <- rowSums(cpm(cts)>10.0) >= 3 #MSP keeps tags which are 10+ in at least 3 libraries
  cts<- cts[keep,]
  dim(cts)

  y <- DGEList(counts=cts ,group=meta$group)
  y <- calcNormFactors(y)

  y<-estimateCommonDisp(y)
  # TagwiseDisp n-value should be close to: 50/(#samples - #groups) = 50/(36-6) = 50/30 =1.666667
  y <- estimateTagwiseDisp(y, prior.n=16)
  
  group<-levels(as.factor(meta$group))		
  et<-exactTest(y, pair=c(group[1],group[2]))
  tab<-summary(de <- decideTestsDGE(et, p=0.0001, adjust="BH"))
  n<-tab[1]+tab[3]
  detags <- rownames(y)[as.logical(de)]
  
  
  plotSmear(et, de.tags=detags, main="DGE Exact Test")
  abline(h = c(-2, 2), col = "blue")
  abline(h = c(-4, 4), col = "blue")
  
  plotMDS.DGEList(y , main = "MDS Plot for Count Data", col=colors)
  textcol<-colors[seq(1,length(colors),5)]
  legend("bottomright",legend=group,text.col= textcol)
  
  plotBCV(y, main="BCV plot")
  
  meanVarPlot <- plotMeanVar(estimateCommonDisp(y) , show.raw.vars=TRUE,
                             show.tagwise.vars=TRUE,
                             show.binned.common.disp.vars=FALSE,
                             show.ave.raw.vars=FALSE , NBline=TRUE,
                             nbins=100,
                             pch=16,
                             xlab="Mean Expresion (Log10)",
                             ylab="Variance (Log10)",
                             main="Mean-Variance Plot")
  
  #positive FC: higher in group2 than 1 
  return(list(topTags(et, n=n),et,cts))
}


pdf(file=pdffile,width=11, height=8)


#evaluate skin vs slime gland

par(mfrow=c(2,2),oma=c(1,1,2,0))
SKvSG<-runEdgeRpwise(quantfiles, cols,sample_meta)
title("Skin vs Slime Gland",outer=T)
write.table(SKvSG[1], "SKvSG_edgeR_output.txt", sep="\t")

res<-data.frame(SKvSG[[2]]$table)
toptags<-data.frame(SKvSG[[1]])

#label edit
rownames(toptags)<-rownames(toptags) %>% gsub("ind4_t."," ",.) %>%  gsub("stoutii_","S ",.)
rownames(res)<-rownames(res) %>% gsub("ind4_t."," ",.) %>%  gsub("stoutii_","S ",.)
rownames(SKvSG[[3]])<-rownames(SKvSG[[3]]) %>% gsub("ind4_t."," ",.) %>%  gsub("stoutii_","S ",.)

toptmp<-toptags %>% tibble::rownames_to_column(var="tag")

#read trees

#read trees and match tips to assembly tags
#atr<-read.tree("4_4_alpha.tre")
#atr<-read.tree("4_1_C.tre")
#atr<-read.tree("../trees/alpha.txt")
atr<-ape::ladderize(read.tree("../trees/alpha_11_8.tree")) #new, full alpha tree

atr$tip.label<-unlist(lapply(atr$tip.label,function(x) strsplit(x,"..",fixed=T)[[1]][[1]])) %>% gsub("'","",.,fixed=T)
atr$tip.label<-gsub("-.*","",atr$tip.label) 
atr$tip.label<-gsub("(Estoutii_\\d+)_\\d+","\\1",atr$tip.label,perl=T)
#label edit
atr$tip.label<-atr$tip.label %>% gsub("ind4_t."," ",.) %>%  gsub("stoutii_","S ",.)

par(mfrow=c(1,1),oma=c(0,0,0,0))
#plot(atr)

#gtr<-read.tree("4_1_gamma.tre")
#gtr<-read.tree("../trees/gamma.txt")
gtr<-ladderize(read.tree("../trees/gamma_11_8_ed1.fas.tre"))
gtr<-ladderize(phytools::reroot(gtr,node=118))

gtr$tip.label<-unlist(lapply(gtr$tip.label,function(x) strsplit(x,"..",fixed=T)[[1]][[1]])) %>% gsub("'","",.,fixed=T)
gtr$tip.label<-gsub("-.*","",gtr$tip.label) 
gtr$tip.label<-gsub("(Estoutii_\\d+)_\\d+","\\1",gtr$tip.label,perl=T)
#label edit
gtr$tip.label<-gtr$tip.label %>% gsub("ind4_t."," ",.) %>%  gsub("stoutii_","S ",.)

pa<-ggtree(atr,right=T,layout='roundrect') + geom_tiplab(color="aquamarine3",size=1) + xlim(0,6)
pg<-ggtree(gtr,right=T,layout='roundrect') + geom_tiplab(color="darkmagenta",size=1) + xlim(0,6)



my_tagsA<-atr$tip.label[atr$tip.label %in% rownames(res)]
my_tagsG<-gtr$tip.label[gtr$tip.label %in% rownames(res)]



#my_tags<-c("E_goslinei_t.44783","E_goslinei_t.102104") #your genes of interest by name
#my_tags<-rownames(res)[res$logCPM>12] #your genes of interest by abundance
#my_tags<-rownames(res)[res$PValue<0.001] #your genes of interest by significance
highlight_dfA <- res[my_tagsA,]%>% tibble::rownames_to_column(var="tag") %>% left_join(toptmp) %>% mutate(tissue=ifelse(FDR<0.0001,ifelse(logFC>0,"skin","slime gland"),NA)) %>% tibble::column_to_rownames(var="tag")
highlight_dfG <- res[my_tagsG,] %>% tibble::rownames_to_column(var="tag") %>% left_join(toptmp) %>% mutate(tissue=ifelse(FDR<0.0001,ifelse(logFC>0,"skin","slime gland"),NA)) %>% tibble::column_to_rownames(var="tag")
hlA.EGind4<- tibble::rownames_to_column(highlight_dfA,var="taxa")
hlG.EGind4<- tibble::rownames_to_column(highlight_dfG,var="taxa")

#save alpha and gamma rows from edgeR
tmp2<-rbind.data.frame(highlight_dfA,highlight_dfG) %>% tibble::rownames_to_column(var="tag")
out<-dplyr::left_join(tmp2,toptmp)
#add counts for tags of interest
mycts<-data.frame(SKvSG[[3]]) %>% tibble::rownames_to_column(var="tag")
out<-left_join(out,mycts)

write.table(out, "EGind4_alpha_gamma_edgeR.txt", sep="\t")


#prune out tips from tree not in filtered tags
#keep other species for context
if(file.exists("../ES_ind1/ESind1_alpha_gamma_edgeR.txt")){
  tmp<-read.table("../ES_ind1/ESind1_alpha_gamma_edgeR.txt")
  bothtags<-c(out$tag,tmp$tag)
  tophagA<-which(atr$tip.label %in% bothtags); otherhagA<-which(atr$tip.label %in% tmp$tag)
  nothagA<-grep("EG|ES",atr$tip.label,invert=T)
  atrFilt<-ape::keep.tip(atr,c(tophagA,nothagA))
  tophagG<-which(gtr$tip.label %in% bothtags); otherhagG<-which(gtr$tip.label %in% tmp$tag)
  nothagG<-grep("EG|ES",gtr$tip.label,invert=T)
  sort(c(otherhagG, nothagG))
  gtrFilt<-ape::keep.tip(gtr,c(tophagG,nothagG))
  
  #replot
  pa<-ggtree(atrFilt,right=T,layout='roundrect') +  geom_tiplab(aes(subset=grepl("EG",label)==F & isTip),color="black",size=1) + geom_tiplab(aes(subset=grepl("EG",label)),color="aquamarine3",size=1)+ xlim(0,6) #+  geom_point2(aes(subset = as.numeric(sub(".*/", "", label))>90 & !isTip),size=2)
  pg<-ggtree(gtrFilt,right=T,layout='roundrect')  + geom_tiplab(aes(subset=grepl("EG",label)==F & isTip),color="black",size=1) + geom_tiplab(aes(subset=grepl("EG",label)),color="darkmagenta",size=1) + xlim(0,6)# +    geom_point2(aes(subset = as.numeric(sub(".*/", "", label))>90 & !isTip),size=2)
 
  #replot better
   #prune off outgroups, then get new node numbers, use to collapse outgroup clade in ggtree  & show nodes based on support < 75
  atrFilt1<-ape::extract.clade(atrFilt,node=116)
  pa<-ggtree(atrFilt1,right=T,layout='roundrect')%<+% hlA.EGind4 + geom_tippoint(aes(shape=tissue),color="aquamarine3",size=2) + scale_shape_manual(values=c(1,2)[c("skin","slime gland") %in% hlA.EGind4$tissue],na.translate=F) +  
    geom_tiplab(aes(subset=grepl("EG",label)==F & isTip),color="black",size=2,offset=0.1) + 
    geom_tiplab(aes(subset=grepl("EG",label)),color="aquamarine3",size=2,offset=0.1)+ xlim(0,4) +  
    geom_point2(aes(subset = as.numeric(sub(".*/", "", label))> 75 & !isTip),size=1) +
    geom_point2(aes(subset = as.numeric(sub("/.*", "", label))> 75 & !isTip),size=1)
  pa<-scaleClade(pa,94, 0.1) %>% collapse(94, 'min',fill="gray") +theme(legend.position = "none")
  
  # repeat with gamma
  gtrFilt1<-ape::extract.clade(gtrFilt,node=111)
  pg<-ggtree(gtrFilt1,right=T,layout='roundrect')%<+% hlG.EGind4 + geom_tippoint(aes(shape=tissue),color="darkmagenta",size=2) + scale_shape_manual(values=c(1,2)[c("skin","slime gland") %in% hlG.EGind4$tissue],na.translate=F) +  
    geom_tiplab(aes(subset=grepl("EG",label)==F & isTip),color="black",size=2,offset=0.1) + 
    geom_tiplab(aes(subset=grepl("EG",label)),color="darkmagenta",size=2,offset=0.1)+ xlim(0,4) +  
    geom_point2(aes(subset = as.numeric(sub(".*/", "", label))> 75 & !isTip),size=1) +
    geom_point2(aes(subset = as.numeric(sub("/.*", "", label))> 75 & !isTip),size=1)
  pg<-pg %>% rotate(78) %>% rotate(79)%>% scaleClade(112, 0.1) %>% collapse(112, 'min',fill="gray") +theme(legend.position = "none")
  
  }


highlight_df2 <- res[rownames(res) %in% rownames(toptags),]
highlight_nonsig <- res[!(rownames(res) %in% rownames(toptags)),]
plotMA_EG<-ggplot(res,aes(x=logCPM,y=logFC)) + 
  labs(x="Mean transcript abundance (log CPM)", y="Differential expression (log FoldChange)") +
  geom_point(data=highlight_nonsig, aes(x=logCPM,y=logFC),col="black",alpha=0.2,size=2) +
  geom_point(data=highlight_df2,aes(x=logCPM,y=logFC,col="red"),alpha=0.2,size=2)  + 
  geom_point(data=highlight_dfA, aes(x=logCPM,y=logFC,shape=tissue),size=3,color="aquamarine3") + 
  geom_point(data=highlight_dfG, aes(x=logCPM,y=logFC,shape=tissue),size=3,color="darkmagenta")  +
  ggrepel::geom_text_repel(data=highlight_dfA,aes(x=logCPM,y=logFC,label=rownames(highlight_dfA)),size = 2) +
  ggrepel::geom_text_repel(data=highlight_dfG,aes(x=logCPM,y=logFC,label=rownames(highlight_dfG)),size = 2) + theme_classic() +theme(legend.position = "none") +scale_shape_manual(values=c(1,2))

rightplots<-plot_grid(pa,pg,ncol=1,align="v",labels=c("B","C"))
cp<-plot_grid(plotMA_EG,rightplots,ncol=2,labels=c("A",""),rel_widths = c(1.5,1))
ggdraw(add_sub(cp,"Figure. Differentially expressed transcripts for alpha and gamma in EGind4. Transcript abundance and differential expression are shown, \nwith transcripts with significant expression differences in red (FDR < 0.0001). Negative logFC values indicate higher expression in slime gland compared to skin. \nTranscripts corresponding to alpha and gamma are annotated. Phylogenetic identification of alpha (B) and gamma (C). Node support > 75% indicated by black points.",x=0,hjust=0,size=10))
dev.off()

MA_EG<-cowplot::as_grob(plotMA_EG)
save(atrFilt1,gtrFilt1,hlA.EGind4,hlG.EGind4,MA_EG,file="EGind4.RData")

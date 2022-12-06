library(ape)
library(ggtree)
library(ggplot2); library(cowplot)
library(tidyr); library(dplyr)
library(this.path)
dir <-this.dir()
setwd(dir)

egtags<-read.table("../EG_ind4/EGind4_alpha_gamma_edgeR.txt") %>% select(tag,tissue)
estags<-read.table("../ES_ind1/ESind1_alpha_gamma_edgeR.txt") %>% select(tag,tissue)
hagtags<-rbind.data.frame(egtags,estags)

atr<-ape::ladderize(read.tree("alpha_11_8.tree")) #new, full alpha tree
atr<-ladderize(phytools::reroot(atr,node=123,position = 0.1))
atr<-extract.clade(atr, node=137)
atr$tip.label<-unlist(lapply(atr$tip.label,function(x) strsplit(x,"..",fixed=T)[[1]][[1]])) %>% gsub("'","",.,fixed=T)
atr$tip.label<-gsub("-.*","",atr$tip.label) 
atr$tip.label<-gsub("(Estoutii_\\d+)_\\d+","\\1",atr$tip.label,perl=T)
#label edit
atr$tip.label<-atr$tip.label %>% gsub("ind4_t."," ",.) %>%  gsub("stoutii_","S ",.)

pdf(file="supp_trees2.pdf",height=11,width=8,useDingbats = F)
pa<-ggtree(atr,right=T,layout='rect')%<+% hagtags + 
  geom_tippoint(aes(shape=tissue),color="aquamarine3",size=3) + 
  scale_shape_manual(values=c(1,2)[c("skin","slime gland") %in% hagtags$tissue],na.translate=F) +  
  geom_tiplab(offset=0.1,size=2)+ xlim(0,6) +  
  geom_point2(aes(subset = as.numeric(sub(".*/", "", label))> 70 & !isTip),size=0.5,color="gray") +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))> 70 & !isTip),size=0.5,color="gray") +   
  geom_point2(aes(subset = as.numeric(sub(".*/", "", label))> 80 & !isTip),size=1,pch=1) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))> 80 & !isTip),size=1,pch=1) +
  geom_point2(aes(subset = as.numeric(sub(".*/", "", label))> 90 & !isTip),size=2,pch=1) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))> 90 & !isTip),size=2,pch=1) +
  theme(legend.position = "top") +labs(shape="Alpha expression enriched:")
#+geom_tiplab(aes(subset=grepl("EGind",label)==F & grepl("Estout",label)==F & isTip),color="black",size=2,offset=0.1) +geom_tiplab(aes(subset=grepl("EGind",label)),color="aquamarine3",size=2,offset=0.1)+ geom_tiplab(aes(subset=grepl("Estout",label)),color="aquamarine3",size=2,offset=0.1)
#pa
ggdraw(add_sub(pa, label="Figure. Alpha phylogeny. Transcripts enriched in skin or slime gland indicated by circles or triangles, respectively. \nNode symbol size represents supports of (70, 80, 90)% bootstrap and/or aLRT.",x=0,hjust=0,size=10))
#dev.off()

library(phangorn)


gtr<-ladderize(read.tree("gamma_11_8_ed1.fas.tre"))
gtr<-ladderize(phytools::reroot(gtr,node=129,position = 0.05))
gtr<-drop.tip(gtr,tip=unlist(Descendants(gtr,node=190,type="tips")))
gtr$tip.label<-unlist(lapply(gtr$tip.label,function(x) strsplit(x,"..",fixed=T)[[1]][[1]])) %>% gsub("'","",.,fixed=T)
gtr$tip.label<-gsub("-.*","",gtr$tip.label) 
gtr$tip.label<-gsub("(Estoutii_\\d+)_\\d+","\\1",gtr$tip.label,perl=T)
#label edit
gtr$tip.label<-gtr$tip.label %>% gsub("ind4_t."," ",.) %>%  gsub("stoutii_","S ",.)

pg<-ggtree(gtr,right=T,layout='rect')%<+% hagtags + 
  geom_tippoint(aes(shape=tissue),color="darkmagenta",size=3) + 
  scale_shape_manual(values=c(1,2)[c("skin","slime gland") %in% hagtags$tissue],na.translate=F) +  
  geom_tiplab(offset=0.1,size=2)+ xlim(0,6) +  
  geom_point2(aes(subset = as.numeric(sub(".*/", "", label))> 70 & !isTip),size=0.5,color="gray") +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))> 70 & !isTip),size=0.5,color="gray") +   
  geom_point2(aes(subset = as.numeric(sub(".*/", "", label))> 80 & !isTip),size=1,pch=1) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))> 80 & !isTip),size=1,pch=1) +
  geom_point2(aes(subset = as.numeric(sub(".*/", "", label))> 90 & !isTip),size=2,pch=1) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))> 90 & !isTip),size=2,pch=1) +
  theme(legend.position = "top") +labs(shape="Gamma expression enriched:")
#+geom_tiplab(aes(subset=grepl("EGind",label)==F & grepl("Estout",label)==F & isTip),color="black",size=2,offset=0.1) +geom_tiplab(aes(subset=grepl("EGind",label)),color="aquamarine3",size=2,offset=0.1)+ geom_tiplab(aes(subset=grepl("Estout",label)),color="aquamarine3",size=2,offset=0.1)
#pa
ggdraw(add_sub(pg, label="Figure. Gamma phylogeny. Transcripts enriched in skin or slime gland indicated by circles or triangles, respectively. \nNode symbol size represents supports of (70, 80, 90)% bootstrap and/or aLRT.",x=0,hjust=0,size=10))
#.
#dev.off()

#make composite figure with MA plots and trees for both species
load("../ES_ind1/ESind1.RData")
load("../EG_ind4/EGind4.RData")

hlA<-rbind.data.frame(hlA.EGind4,hlA.ESind1)
  pa<-ggtree(atrFilt1,right=T,layout='roundrect')%<+% hlA + geom_tippoint(aes(shape=tissue),color="aquamarine3",size=2) + scale_shape_manual(values=c(1,2)[c("skin","slime gland") %in% hlA$tissue],na.translate=F) +  
  geom_tiplab(aes(subset=grepl("E[G|s]",label)==F & isTip),color="black",size=2,offset=0.1) + 
  geom_tiplab(aes(subset=grepl("E[G|s]",label)),color="aquamarine3",size=2,offset=0.1)+ xlim(0,4) +  
  geom_point2(aes(subset = as.numeric(sub(".*/", "", label))> 75 & !isTip),size=1) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))> 75 & !isTip),size=1)
pa<-scaleClade(pa,94, 0.1) %>% collapse(94, 'min',fill="gray") +theme(legend.position = "none")

# repeat with gamma
hlG<-rbind.data.frame(hlG.EGind4,hlG.ESind1)

pg<-ggtree(gtrFilt1,right=T,layout='roundrect')%<+% hlG + geom_tippoint(aes(shape=tissue),color="darkmagenta",size=2) + scale_shape_manual(values=c(1,2)[c("skin","slime gland") %in% hlG$tissue],na.translate=F) +  
  geom_tiplab(aes(subset=grepl("E[G|s]",label)==F & isTip),color="black",size=2,offset=0.1) + 
  geom_tiplab(aes(subset=grepl("E[G|s]",label)),color="darkmagenta",size=2,offset=0.1)+ xlim(0,4) +  
  geom_point2(aes(subset = as.numeric(sub(".*/", "", label))> 75 & !isTip),size=1) +
  geom_point2(aes(subset = as.numeric(sub("/.*", "", label))> 75 & !isTip),size=1)
pg<-pg %>% rotate(78) %>% rotate(79)%>% scaleClade(112, 0.1) %>% collapse(112, 'min',fill="gray") +theme(legend.position = "none")

topplots<-plot_grid(MA_EG,MA_ES,align="h",labels=c("A","B"),ncol=2)
botplots<-plot_grid(pa,pg,ncol=2,align="h",labels=c("C","D"))
cp<-plot_grid(topplots,botplots,ncol=1,rel_heights = c(1,1))
ggdraw(add_sub(cp,"Figure. Differentially expressed transcripts for alpha and gamma. Transcript abundance and differential expression are shown for\n EGind4 (A) and ESind1 (B), with transcripts with significant expression differences in red (FDR < 0.0001). Negative logFC values \nindicate higher expression in slime gland compared to skin. Transcripts corresponding to alpha and gamma are annotated. \nPhylogenetic placement of alpha (C) and gamma (D) transcripts. Node support > 75% indicated by black points.",x=0,hjust=0,size=10))
dev.off()
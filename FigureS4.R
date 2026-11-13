rm(list=ls())
setwd("") # Set this to the local 05_figures folder before running
library(ggplot2)
library(dplyr)
library(svglite)
library(tidytext)
library(cowplot)

dir.create("output",showWarnings=FALSE)

####################Figure S4A: GTEx tissue expression####################
# 1. Draw the 115-gene tissue profile from previously calculated sampling summaries
x <- read.csv("data/GTEx_sampling.csv",stringsAsFactors=FALSE)
stopifnot(all(c("GeneSet","Tissue","Group","Mean","Lower","Upper") %in% names(x)))
x <- x %>% filter(GeneSet=="115")
if(nrow(x)==0) stop("GTEx 115-gene plotting data missing: data/GTEx_sampling.csv")
x$Tissue <- factor(x$Tissue,levels=unique(x$Tissue))
x$Group <- factor(x$Group,levels=c("Control","Hotspot","NDD_associated"))
stopifnot(!anyNA(x$Group),all(x$Mean>=0 & x$Mean<=100))
colors <- c("Control"="#377EB8","Hotspot"="#E69F00","NDD_associated"="#D73027")
pS4a <- ggplot(x,aes(Tissue,Mean,color=Group,group=Group))+
  geom_ribbon(data=subset(x,Group!="Hotspot"),aes(ymin=Lower,ymax=Upper,fill=Group),
              alpha=0.13,color=NA)+geom_line(linewidth=0.4)+geom_point(size=0.7)+
  scale_color_manual(values=colors)+scale_fill_manual(values=colors)+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))+
  labs(x=NULL,y="Proportion of expressed genes (%)",color=NULL,fill=NULL)+
  theme_classic(base_size=10)+
  theme(axis.text.x=element_text(size=6,angle=45,hjust=1),legend.position="top")
print(pS4a)
ggsave("output/FigureS4A_GTEx.svg",pS4a,width=22,height=15,units="cm")

####################Figure S4B: scRNA-seq class-specific module scores####################
# 1. Read 21 class-by-cell-type tests; FDR was computed across all 21 comparisons
s <- read.csv("data/scRNA_bubble.csv",stringsAsFactors=FALSE)
stopifnot(all(c("Class","Celltype","Mean_Z","P","FDR") %in% names(s)),nrow(s)==21)
stopifnot(all(is.finite(s$Mean_Z)),all(s$FDR>=0 & s$FDR<=1))
class_order <- c("Both","Denovonear-only","3D-CLUMP-only")
cell_order <- c("Excitatory neuron","Radial glia","Inhibitory neuron",
                "Fibroblast","Microglia","Vascular-associated cell","OPC-like cells")
s$Class <- factor(s$Class,levels=rev(class_order))
s$Celltype <- factor(s$Celltype,levels=cell_order)
stopifnot(!anyNA(s$Class),!anyNA(s$Celltype),
          !anyDuplicated(s[c("Class","Celltype")]))
s$Significance <- ifelse(s$FDR<0.001,"***",ifelse(s$FDR<0.01,"**",ifelse(s$FDR<0.05,"*","")))
s$Bubble_size <- pmin(-log10(pmax(s$FDR,.Machine$double.xmin)),10)
z_limit <- max(0.5,max(abs(s$Mean_Z)))
pS4b <- ggplot(s,aes(x=Celltype,y=Class))+
  geom_point(aes(size=Bubble_size,fill=Mean_Z),shape=21,color="grey25",stroke=0.35)+
  geom_text(aes(label=Significance),color="black",size=3.3,fontface="bold")+
  scale_fill_gradient2(low="#4575B4",mid="#F7F7F7",high="#D73027",
                       midpoint=0,limits=c(-z_limit,z_limit),name="Mean Z-score")+
  scale_size_continuous(range=c(2.5,12),limits=c(0,10),breaks=c(1,2,3,5,10),
                        name=expression(-log[10]("BH-adjusted P")))+
  labs(x=NULL,y=NULL)+theme_classic(base_size=11)+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank(),
        axis.text.x=element_text(angle=35,hjust=1,size=9),legend.position="right")
print(pS4b)
ggsave("output/FigureS4B_scRNA_bubble.svg",pS4b,width=16,height=8,units="cm")

####################Supplementary GO, KEGG and DisGeNET enrichment####################
# 1. Plot selected significant terms from the original analyses; no API access needed
z <- read.csv("data/Enrichment_plot.csv",stringsAsFactors=FALSE)
required <- c("Analysis","Database","Class","Description","GeneRatio_num","Count","FDR")
stopifnot(all(required %in% names(z)),nrow(z)>0)
stopifnot(all(z$GeneRatio_num>=0),all(z$Count>=0),all(z$FDR>0 & z$FDR<=1))
# 2. Produce 1x3 dot-plot sets for each analysis in the input table
for(analysis in unique(z$Analysis)){
  d <- z %>% filter(Analysis==analysis)
  if(analysis=="methods") db_order <- c("GO_BP","KEGG","DisGeNET")
  else if(analysis=="annotation") db_order <- c("GO_MF","KEGG","DisGeNET")
  else stop("Unknown enrichment analysis: ",analysis)
  plots <- lapply(db_order,function(db){
    g <- d %>% filter(Database==db)
    if(nrow(g)==0) return(ggplot()+theme_void()+labs(title=db))
    g$Class <- factor(g$Class,levels=unique(g$Class))
    g$Term <- tidytext::reorder_within(g$Description,g$GeneRatio_num,g$Class)
    ggplot(g,aes(x=GeneRatio_num,y=Term))+
      geom_point(aes(size=Count,fill=FDR),shape=21,color="black",stroke=0.3)+
      scale_size(range=c(3,8),name="Count")+
      scale_fill_gradient(low="#2166AC",high="#B2182B",name="BH-adjusted P")+
      tidytext::scale_y_reordered()+
      facet_grid(Class~.,scales="free_y",space="free_y")+
      labs(x="Gene ratio",y=NULL,title=db)+theme_classic(base_size=10)+
      theme(strip.text=element_text(face="bold"),plot.title=element_text(face="bold"),
            axis.text.y=element_text(size=8))
  })
  p <- cowplot::plot_grid(plotlist=plots,nrow=1,align="h",axis="tb")
  print(p)
  ggsave(paste0("output/FigureS4_enrichment_",analysis,".svg"),p,
         width=52,height=18,units="cm")
}

rm(list=ls())
setwd("") # Set to this repository's 04_figures directory before Source
####################Figure S6C: algorithm-specific continuous outcomes####################
# 1. Packages and data
library(ggplot2)
library(svglite)
DATA <- "data"
OUT <- "output"
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
f <- file.path(DATA,"FigureS6C_HC3.csv")
if(!file.exists(f))stop("Missing: ",f)
d <- read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
required <- c("phenotype","label","beta","p_raw","p_adj","status","n_ThreeD","n_Near")
if(length(setdiff(required,names(d))))stop("FigureS6C_HC3.csv missing required columns")
SEL <- c("scq_total","rbsr_total","rbsr_compulsive","rbsr_stereotyped","age_onset_mos")
if(nrow(d)!=5L||anyDuplicated(d$phenotype)||!setequal(d$phenotype,SEL))stop("Expected five S6C traits")
d <- d[match(SEL,d$phenotype),,drop=FALSE]
# 2. Display ORIGINAL raw HC3 P values, not the BH5 values.
p_label <- function(x)ifelse(is.finite(x),ifelse(x<.001,"p < 0.001",sprintf("p = %.3f",x)),"NT")
d$label <- factor(ifelse(is.na(d$label)|!nzchar(d$label),d$phenotype,d$label),levels=rev(ifelse(is.na(d$label)|!nzchar(d$label),d$phenotype,d$label)))
d$cell <- p_label(d$p_raw)
lim <- suppressWarnings(max(abs(d$beta),na.rm=TRUE));if(!is.finite(lim)||lim<.2)lim<-.2;lim<-ceiling(lim*10)/10
p <- ggplot(d,aes(x="3D-clump vs Denovonear",y=label))+
  geom_tile(aes(fill=beta),colour="white",linewidth=.8,width=.95,height=.85)+
  geom_text(aes(label=cell),size=3.1,colour="#202020")+
  scale_fill_gradient2(low="#4C83A8",mid="white",high="#B23F36",midpoint=0,limits=c(-lim,lim),na.value="grey95",name="Standardized β")+
  scale_x_discrete(position="top")+
  labs(title="(C) Continuous phenotypes: OLS + HC3",x=NULL,y=NULL)+
  theme_minimal(base_size=10)+theme(panel.grid=element_blank(),axis.text.x=element_text(colour="black",size=8),axis.text.y=element_text(colour="black",size=9),legend.position="right")
print(p)
ggsave(file.path(OUT,"FigureS6C.svg"),plot=p,width=4,height=3,device=svglite::svglite)
ggsave(file.path(OUT,"FigureS6C.pdf"),plot=p,width=4,height=3,device="pdf")
message("Figure S6C completed. Source BH correction: five selected traits; plotted labels: raw P.")

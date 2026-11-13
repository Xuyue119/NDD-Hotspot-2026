library(ggplot2)
library(svglite)
####################Figure 2B-E####################
dir.create("output",showWarnings=FALSE)
red <- "#BF1D2D"; blue <- "#3D549D"
star <- function(p) ifelse(is.na(p),"",ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",""))))
save_panel <- function(p,name,w=8.5,h=7.5){print(p);ggsave(file.path("output",paste0(name,".svg")),p,width=w,height=h,units="cm");ggsave(file.path("output",paste0(name,".pdf")),p,width=w,height=h,units="cm")}
base_theme <- theme_classic(base_size=10)+theme(axis.text=element_text(color="black"),axis.title=element_text(color="black"),axis.line=element_line(linewidth=0.5),axis.ticks=element_line(linewidth=0.5))
####################1. Figure 2B-C: burden forest plots####################
burden <- read.csv("data/Figure2_BC.csv",check.names=FALSE)
plot_forest <- function(panel,levels,xlim,xbreaks){
  d <- subset(burden,Panel==panel)
  if(nrow(d)!=3 || !setequal(d$Category,levels)) stop(paste("Unexpected burden categories in panel",panel))
  d$Category <- factor(d$Category,levels=levels)
  p <- ggplot(d,aes(x=RR,y=Category))+geom_vline(xintercept=1,linetype=2,color="#9A9A9A",linewidth=0.45)+
    geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=0.1,linewidth=0.6,color="#333333")+
    geom_point(size=2.1,color="#333333")+
    geom_text(data=subset(d,!is.na(p_adj_BH) & p_adj_BH<0.05),aes(label=star(p_adj_BH)),nudge_y=0.24,size=3.1)+
    scale_x_continuous(breaks=xbreaks,expand=c(0,0))+coord_cartesian(xlim=xlim)+
    labs(x="Rate ratio (Proband/Sibling)",y=NULL)+base_theme+
    theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
  save_panel(p,paste0("Figure2",panel),6,5)
}
plot_forest("B",c("MIS","LGD","SYN"),c(0.90,1.51),seq(0.9,1.5,0.1))
plot_forest("C",c("CADD ≥ 30","20 ≤ CADD < 30","CADD < 20"),c(0.75,1.33),seq(0.8,1.3,0.1))
####################2. Shared score ECDF input and annotations####################
ecdf <- read.csv("data/FigureS2_Predictor_ECDF.csv",check.names=FALSE)
stats <- read.csv("data/FigureS2_Plot_statistics.csv",check.names=FALSE)
ecdf$Group <- factor(ecdf$Group,levels=c("Proband","Sibling"))
panel_star <- function(panel){s <- stats$p_adj_BH[stats$Panel==panel];if(length(s)!=1) stop(paste("Missing/duplicated adjusted P for",panel));star(s)}
####################3. Figure 2D: CADD ECDF####################
d <- subset(ecdf,Predictor=="CADD")
if(!nrow(d)) stop("No CADD ECDF data")
p2d <- ggplot(d,aes(x=Score,y=Cumulative_proportion,color=Group,group=Group))+
  geom_step(linewidth=0.75,direction="hv")+
  geom_vline(xintercept=c(20,30),linetype=2,linewidth=0.45,color="#8A8A8A")+
  annotate("text",x=4,y=0.94,label=panel_star("Figure2D"),size=4)+
  scale_color_manual(values=c(Proband=red,Sibling=blue))+
  scale_x_continuous(breaks=seq(0,50,10),expand=c(0,0))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.25),expand=c(0,0))+
  labs(x="CADD PHRED score",y="Cumulative proportion",color=NULL)+base_theme+
  theme(legend.position="top",legend.justification="left")
save_panel(p2d,"Figure2D")
####################4. Figure 2E: gene-level CADD by dnMIS recurrence####################
gene <- read.csv("data/Figure2E_gene_CADD.csv",check.names=FALSE)
gene$Category <- factor(gene$Category,levels=c("1","2","3","4+"))
if(!nrow(gene) || anyNA(gene$Category)) stop("Invalid gene recurrence input")
comparisons <- data.frame(Panel="Figure2E",Comparison=c("1 vs 2","1 vs 3","1 vs 4+"),x1=1,x2=c(2,3,4),y=c(41,43.5,46))
sig <- merge(comparisons,stats[,c("Panel","Comparison","p_adj_BH")],by=c("Panel","Comparison"),all.x=TRUE,sort=FALSE)
if(anyNA(sig$p_adj_BH)) stop("Missing Figure 2E adjusted P values")
sig$label <- star(sig$p_adj_BH); sig <- subset(sig,label!="")
p2e <- ggplot(gene,aes(x=Category,y=CADD_median,fill=Category))+
  geom_violin(width=0.85,trim=FALSE,color="#595959",linewidth=0.45)+
  geom_boxplot(width=0.14,outlier.shape=NA,fill="white",color="black",linewidth=0.45)+
  geom_segment(data=sig,aes(x=x1,xend=x2,y=y,yend=y),inherit.aes=FALSE,linewidth=0.45)+
  geom_text(data=sig,aes(x=(x1+x2)/2,y=y+0.8,label=label),inherit.aes=FALSE,size=4)+
  scale_fill_manual(values=c("1"="#D9D9D9","2"="#BDBDBD","3"="#969696","4+"="#737373"))+
  scale_y_continuous(limits=c(0,48),breaks=seq(0,40,10),expand=c(0,0))+
  labs(x="Number of dnMIS variants per gene",y="Median CADD PHRED score per gene")+
  base_theme+theme(legend.position="none")
save_panel(p2e,"Figure2E")

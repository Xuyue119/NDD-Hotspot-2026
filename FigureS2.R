library(ggplot2)
library(svglite)
####################Figure S2A-D####################
dir.create("output",showWarnings=FALSE)
red <- "#BF1D2D";blue <- "#3D549D"
star <- function(p) ifelse(is.na(p),"",ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",""))))
save_panel <- function(p,name){print(p);ggsave(file.path("output",paste0(name,".svg")),p,width=8.5,height=7.5,units="cm");ggsave(file.path("output",paste0(name,".pdf")),p,width=8.5,height=7.5,units="cm")}
base_theme <- theme_classic(base_size=10)+theme(axis.text=element_text(color="black"),axis.title=element_text(color="black"),axis.line=element_line(linewidth=0.5),axis.ticks=element_line(linewidth=0.5))
####################1. Figure S2A: matched recurrence####################
d <- read.csv("data/FigureS2A_recurrence.csv",check.names=FALSE)
d <- d[match(c("1","2","3","4+"),d$Category),]
if(nrow(d)!=4 || anyNA(d$Category)) stop("Expected recurrence categories 1, 2, 3, 4+")
d$x <- 1:4
p <- ggplot()+
  geom_errorbar(data=d,aes(x=x-0.07,ymin=Proband_CI_low,ymax=Proband_CI_high),width=0.08,linewidth=0.6,color=red)+
  geom_point(data=d,aes(x=x-0.07,y=Proband_median),shape=16,size=1.8,color=red)+
  geom_point(data=d,aes(x=x+0.07,y=Sibling_gene_count),shape=15,size=1.8,color=blue)+
  annotate("point",x=c(1.55,2.45),y=4630,shape=c(16,15),color=c(red,blue),size=1.8)+
  annotate("text",x=c(1.65,2.55),y=4630,label=c("Proband","Sibling"),hjust=0,size=3.1)+
  scale_x_continuous(breaks=1:4,labels=c("1","2","3","4+"),limits=c(0.6,4.4),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,4000,1000),limits=c(0,4800),expand=c(0,0))+
  labs(x="Number of dnMIS variants per gene",y="Number of genes")+base_theme
save_panel(p,"FigureS2A")
####################2. Figure S2B-D: predictor ECDFs####################
ecdf <- read.csv("data/FigureS2_Predictor_ECDF.csv",check.names=FALSE)
stats <- read.csv("data/FigureS2_Plot_statistics.csv",check.names=FALSE)
ecdf$Group <- factor(ecdf$Group,levels=c("Proband","Sibling"))
plot_ecdf <- function(predictor,panel,xlab,show_legend){
  d <- subset(ecdf,Predictor==predictor)
  p_adj <- stats$p_adj_BH[stats$Panel==paste0("FigureS2",panel)]
  if(!nrow(d) || length(p_adj)!=1) stop(paste("Missing ECDF or adjusted P:",predictor))
  p <- ggplot(d,aes(x=Score,y=Cumulative_proportion,color=Group,group=Group))+
    geom_step(linewidth=0.75,direction="hv")+
    annotate("text",x=0.08,y=0.94,label=star(p_adj),size=4)+
    scale_color_manual(values=c(Proband=red,Sibling=blue))+
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2),expand=c(0,0))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.25),expand=c(0,0))+
    labs(x=xlab,y="Cumulative proportion",color=NULL)+base_theme+
    theme(legend.position=if(show_legend) "top" else "none",legend.justification="left")
  save_panel(p,paste0("FigureS2",panel))
}
plot_ecdf("REVEL","B","REVEL rank score",TRUE)
plot_ecdf("MPC","C","MPC rank score",FALSE)
plot_ecdf("AlphaMissense","D","AlphaMissense rank score",FALSE)

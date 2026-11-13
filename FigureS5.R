rm(list=ls())
setwd("") # Set to 05_figures
####################Figure S5B-D: Both-clustered versus Non-clustered####################
# 1. Packages and inputs
library(ggplot2)
dir.create("output",showWarnings=FALSE)
read_input<-function(f,cols){p<-file.path("data",f);if(!file.exists(p))stop("Missing ",p);d<-read.csv(p,check.names=FALSE,stringsAsFactors=FALSE);m<-setdiff(cols,names(d));if(length(m))stop(f," missing: ",paste(m,collapse=", "));d}
fmtp<-function(p)ifelse(p<.001,"<0.001",sprintf("%.3f",p))
cols<-c(Both="#c75a35",Non="#4C83A8")
forest<-read_input("FigureS5B_forest.csv",c("Score","Label","Beta","CI_low","CI_high","P_BH"))
oe<-read_input("FigureS5CD_oe.csv",c("Figure","Metric","Group","N","Hits","Observed","Expected","Fisher_P","OE_P"))
if(nrow(forest)!=4L||!setequal(oe$Group,names(cols))||!setequal(oe$Figure,c("S5C","S5D")))stop("Figure S5 input mismatch")
# 2. Figure S5B
order_score<-c("CADD_phred","REVEL_rankscore","MPC_rankscore","AlphaMissense_rankscore")
forest<-forest[match(order_score,forest$Score),,drop=FALSE]
if(anyNA(forest$Score))stop("Missing predictor score")
forest$Y<-rev(seq_len(nrow(forest)))
forest$Stars<-ifelse(forest$P_BH<.001,"***",ifelse(forest$P_BH<.01,"**",ifelse(forest$P_BH<.05,"*","")))
p<-ggplot(forest,aes(x=Beta,y=Y))+geom_vline(xintercept=0,linetype="dashed",color="grey65")+geom_segment(aes(x=CI_low,xend=CI_high,yend=Y),linewidth=.7,color="#c75a35")+geom_point(size=2.6,color="#c75a35")+geom_text(aes(label=Stars),nudge_y=.22,size=4)+scale_y_continuous(breaks=forest$Y,labels=forest$Label)+labs(x="Standardized score difference: Both - Non (95% CI)",y=NULL)+theme_classic(base_size=10)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
print(p);ggsave("output/FigureS5B.svg",p,width=9,height=5,units="cm")
# 3. Figure S5C-D
plot_bar<-function(fig,metric,ymax,breaks){d<-oe[oe$Figure==fig & oe$Metric==metric,,drop=FALSE];if(nrow(d)!=2L||!setequal(d$Group,names(cols)))stop("Missing groups for ",fig);d$Group<-factor(d$Group,levels=c("Both","Non"));p<-ggplot(d,aes(x=Group,y=Observed,fill=Group))+geom_col(width=.6)+geom_point(aes(y=Expected),shape=21,fill="white",color="black",size=3.3,stroke=.9)+geom_text(aes(label=sprintf("%d/%d\n%.1f%%",Hits,N,Observed)),vjust=-.4,size=3)+scale_fill_manual(values=cols)+scale_x_discrete(labels=c(Both="Both-clustered",Non="Non-clustered"))+scale_y_continuous(limits=c(0,ymax),breaks=breaks,labels=function(x)paste0(x,"%"),expand=c(0,0))+labs(subtitle=paste0("Fisher p=",fmtp(d$Fisher_P[1]),"\nO/E p=",fmtp(d$OE_P[1])),x=NULL,y=if(fig=="S5C")"CpG proportion of unique dnMIS variants (%)" else "UniProtKB active/binding-site overlap (%)")+theme_classic(base_size=10)+theme(legend.position="none",plot.subtitle=element_text(hjust=.5));print(p);ggsave(paste0("output/Figure",fig,".svg"),p,width=8,height=9,units="cm")}
plot_bar("S5C","CpG",45,seq(0,45,5))
plot_bar("S5D","UniProtKB active/binding",6,seq(0,6,2))

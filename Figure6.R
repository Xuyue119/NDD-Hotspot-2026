rm(list=ls())
setwd("") # Set to this repository's 04_figures directory before Source
####################Figure 6B/C/D: load plotting inputs####################
# 1. Packages; all models and resampling were completed in the original analysis.
library(ggplot2)
library(scales)
library(svglite)
DATA <- "data"
OUT <- "output"
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
read_input <- function(name,required){path<-file.path(DATA,name);if(!file.exists(path))stop("Missing: ",path);x<-read.csv(path,check.names=FALSE,stringsAsFactors=FALSE);m<-setdiff(required,names(x));if(length(m))stop(name," missing columns: ",paste(m,collapse=", "));x}
p_label <- function(x) ifelse(is.finite(x),ifelse(x<.001,"p < 0.001",sprintf("p = %.3f",x)),"NT")
save_plot <- function(p,name,w,h){print(p);ggsave(file.path(OUT,paste0(name,".svg")),plot=p,width=w,height=h,device=svglite::svglite);ggsave(file.path(OUT,paste0(name,".pdf")),plot=p,width=w,height=h,device="pdf")}
GROUPS <- c("Both","Non","Syn")
DISPLAY <- c(Both="Both-clustered variants",Non="Non-clustered variants",Syn="dnSYN variants")
PAL <- c(Both="#B23F36",Non="#4C83A8",Syn="#474554")

####################Figure 6B: HC3 continuous regression heatmap####################
# 1. The plotted raw p-values and standardized betas are from the primary Both vs Non HC3 comparisons.
b <- read_input("Figure6B_HC3.csv",c("phenotype","label","panel","position","beta","p_raw","p_adj","status"))
left <- c("scq_total","rbsr_total","rbsr_compulsive","rbsr_stereotyped","age_onset_mos")
right <- c("dcdq_total","dcdq_fine_motor","dcdq_coordination","vineland_abc","vineland_communication")
if(nrow(b)!=10L||!setequal(b$phenotype,c(left,right))||anyDuplicated(b$phenotype))stop("Incorrect Figure6B phenotype set")
if(any(b$panel[match(left,b$phenotype)]!="Social and repetitive behaviors")||any(b$panel[match(right,b$phenotype)]!="Motor and adaptive functioning"))stop("Figure6B panel mismatch")
b$position <- as.integer(b$position)
if(!identical(b$position[match(c(left,right),b$phenotype)],rep(1:5,2)))stop("Figure6B position mismatch")
b$panel <- factor(b$panel,levels=c("Social and repetitive behaviors","Motor and adaptive functioning"))
b$p_label <- p_label(b$p_raw)
lim <- max(abs(b$beta),na.rm=TRUE);if(!is.finite(lim)||lim<.2)lim<-.2;lim<-ceiling(lim*10)/10
pB <- ggplot(b,aes(x=position,y=1))+
  geom_tile(data=b[is.finite(b$beta),,drop=FALSE],aes(fill=beta),width=.88,height=.48,colour="white",linewidth=.6)+
  geom_text(aes(label=p_label),size=3.1,colour="#202020")+
  geom_text(aes(y=.56,label=label),size=3,lineheight=.95,colour="#202020")+
  facet_wrap(~panel,ncol=1,strip.position="top")+
  scale_fill_gradient2(low="#4C83A8",mid="white",high="#B23F36",midpoint=0,limits=c(-lim,lim),name="Standardized β")+
  scale_x_continuous(limits=c(.45,5.55),breaks=NULL,expand=expansion(mult=0))+
  scale_y_continuous(limits=c(.35,1.38),breaks=NULL,expand=expansion(mult=0))+
  labs(title="(B) Continuous phenotypes",subtitle="Both-clustered vs Non-clustered | OLS with HC3",x=NULL,y=NULL)+
  theme_minimal(base_size=11)+theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),strip.background=element_blank(),strip.text=element_text(size=11,face="bold",hjust=.5),legend.position="right",panel.spacing.y=grid::unit(.25,"cm"),plot.title=element_text(face="bold",size=12),plot.subtitle=element_text(size=9))
save_plot(pB,"Figure6B_continuous",10,7.5)

####################Figure 6C: cognitive ordinal-category composition####################
# 1. Rank 1-10 counts come from the 3-group phenotype RDS; p is the ordinal omnibus LRT.
cog <- read_input("Figure6C_cognitive.csv",c("group","rank","label","n","n_valid","p_omnibus_raw"))
if(nrow(cog)!=30L||!setequal(unique(cog$group),GROUPS)||!setequal(unique(cog$rank),1:10))stop("Figure6C grid inconsistent")
if(anyDuplicated(cog[,c("group","rank")]))stop("Duplicate cognitive cells")
if(any(cog$n<0)||any(cog$n_valid<=0))stop("Invalid cognitive counts")
for(g in GROUPS){z<-cog[cog$group==g,,drop=FALSE];if(sum(z$n)!=unique(z$n_valid))stop("Cognitive denominator mismatch: ",g)}
p_cog <- unique(cog$p_omnibus_raw);if(length(p_cog)!=1L)stop("Cognitive omnibus P not unique")
labels <- c("≤24","25–39","40–54","55–69","70–79","80–89","90–109","110–119","120–129","≥130")
if(!identical(unique(as.character(cog$label[order(cog$rank)])),labels))stop("Cognitive rank labels mismatch")
cog$group <- factor(cog$group,levels=GROUPS)
cog$category <- factor(cog$label,levels=labels)
cols<-grDevices::colorRampPalette(c("#B23F36","#EACFC9","#A6C5DA","#4C83A8"))(10);names(cols)<-labels
pC <- ggplot(cog,aes(x=group,y=n,fill=category))+
  geom_col(position="fill",width=.71,colour="white",linewidth=.25)+
  scale_fill_manual(values=cols,drop=FALSE)+scale_x_discrete(labels=DISPLAY[GROUPS])+
  scale_y_continuous(breaks=seq(0,1,.25),labels=percent_format(),expand=expansion(mult=c(0,0)))+
  labs(title=paste0("(C) Cognitive score (",p_label(p_cog),")"),x=NULL,y="Proportion",fill="Cognitive score")+
  theme_classic(base_size=10)+theme(legend.position="top",axis.text.x=element_text(size=8),legend.text=element_text(size=8))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
save_plot(pC,"Figure6C_cognitive",4,4.2)

####################Figure 6D: crawling Kaplan-Meier cumulative attainment####################
# 1. Replot original complete-case KM steps; omnibus p is from adjusted 3-group Cox model.
km <- read_input("Figure6D_KM.csv",c("group","month","attained","p_Cox_omnibus_raw"))
if(!setequal(unique(km$group),GROUPS))stop("KM missing group")
p_cox <- unique(km$p_Cox_omnibus_raw);if(length(p_cox)!=1L)stop("Cox omnibus P not unique")
for(g in GROUPS){z<-km[km$group==g,,drop=FALSE];z<-z[order(z$month),,drop=FALSE];if(nrow(z)<2L||anyDuplicated(z$month)||any(diff(z$attained) < -1e-10)||any(z$attained<0|z$attained>1))stop("Invalid KM curve: ",g)}
km$group <- factor(km$group,levels=GROUPS)
km <- km[order(km$group,km$month),,drop=FALSE]
pD <- ggplot(km,aes(x=month,y=attained,colour=group,group=group))+
  geom_step(direction="hv",linewidth=.9)+
  scale_colour_manual(values=PAL[GROUPS],labels=DISPLAY[GROUPS],drop=FALSE)+
  scale_x_continuous(breaks=seq(0,36,6),expand=expansion(mult=c(0,0)))+
  scale_y_continuous(breaks=seq(0,1,.25),labels=percent_format(),expand=expansion(mult=c(0,0)))+
  coord_cartesian(xlim=c(0,36),ylim=c(0,1))+
  labs(title=paste0("(D) Crawling (",p_label(p_cox),")"),x="Age (months)",y="Cumulative attainment",colour=NULL)+
  theme_classic(base_size=10)+theme(legend.position="top",axis.text=element_text(colour="black",size=8))
save_plot(pD,"Figure6D_crawl",4,4.2)
message("Figure 6B, 6C, and 6D plots completed. KM is descriptive; Cox p is covariate-adjusted.")

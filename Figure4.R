rm(list=ls())
setwd("") # Set to the GitHub repository's 05_figures folder
####################Figure 4D: scRNA tSNE module score####################
# 1. Packages
library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)

# 2. Files
dir.create("output",showWarnings=FALSE)

expr_file <- "exprMatrix.tsv.gz"
meta_file <- "meta.tsv"
gene_file <- "data/hotspot115_genes.csv"

if(!file.exists(expr_file))
  download.file("https://cells.ucsc.edu/cortex-dev/exprMatrix.tsv.gz",
                expr_file,mode="wb")
if(!file.exists(meta_file))
  download.file("https://cells.ucsc.edu/cortex-dev/meta.tsv",
                meta_file,mode="wb")
if(!file.exists(gene_file)) stop("Missing: ",gene_file)

# 3. Read UCSC cortex-dev data
mat <- fread(expr_file,check.names=FALSE)
meta <- as.data.frame(fread(meta_file),row.names=1,
                      check.names=FALSE)

genes <- sub(".*\\|","",as.character(mat[[1]]))
if(anyNA(genes)||anyDuplicated(genes)||any(genes==""))
  stop("Missing or duplicated gene symbols")

expr <- as.matrix(mat[,-1,with=FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- genes

if(!setequal(colnames(expr),rownames(meta)))
  stop("Expression matrix and metadata cell IDs disagree")

meta <- meta[colnames(expr),,drop=FALSE]
expr <- Matrix(expr,sparse=TRUE)
rm(mat)
gc()

# 4. Seurat object
so <- CreateSeuratObject(
  counts=expr,project="Cortex",
  meta.data=meta,min.cells=0,min.features=0
)
rm(expr,meta)
gc()

# UCSC input is already normalized; do not run NormalizeData().
so[["pct.RB"]] <- PercentageFeatureSet(so,pattern="^RP[LS]")
hb <- intersect(
  c("HBA1","HBA2","HBB","HBD","HBE1",
    "HBG1","HBG2","HBM","HBQ1","HBZ"),
  rownames(so)
)
so[["pct.HB"]] <- if(length(hb)>0)
  PercentageFeatureSet(so,features=hb) else rep(0,ncol(so))

so <- subset(
  so,subset=nFeature_RNA>200 & nFeature_RNA<10000 &
    pct.RB<7.5 & pct.HB<0.05
)

if(ncol(so)!=4238)
  stop("Expected 4,238 cells after QC; found ",ncol(so))

# 5. Log transform and dimensionality reduction
raw_expr <- GetAssayData(so,assay="RNA",layer="counts")
so <- SetAssayData(
  so,assay="RNA",layer="data",
  new.data=log1p(raw_expr)
)
rm(raw_expr)

set.seed(123)
so <- FindVariableFeatures(
  so,selection.method="dispersion",nfeatures=2000
)
so <- ScaleData(so,features=VariableFeatures(so))
so <- RunPCA(so,features=VariableFeatures(so),npcs=30)
so <- FindNeighbors(so,dims=1:8)
so <- FindClusters(so,resolution=0.1,algorithm=1)
so <- RunTSNE(so,dims=1:8,seed.use=123)

# 6. Hotspot gene module score
hotspot <- read.csv(gene_file,stringsAsFactors=FALSE)$Gene
hotspot <- unique(trimws(hotspot))
hotspot <- hotspot[!is.na(hotspot)&nzchar(hotspot)]

if(length(hotspot)!=115)
  stop("Expected 115 hotspot genes")

gene_alias <- c(
  "TAFA4"="FAM19A4","YARS1"="YARS",
  "GALNT17"="WBSCR17","GFUS"="TSTA3",
  "SEPTIN2"="SEPT2","CEMIP2"="TMEM2"
)
idx <- hotspot %in% names(gene_alias)
hotspot[idx] <- unname(gene_alias[hotspot[idx]])
hotspot <- intersect(unique(hotspot),rownames(so))

if(length(hotspot)<5)
  stop("Too few hotspot genes detected")

so <- AddModuleScore(
  so,features=list(hotspot),
  name="hotspot_module_score",
  assay="RNA",slot="data",
  ctrl=100,seed=123
)

# 7. Plot
coords <- Embeddings(so,reduction="tsne")
score <- so@meta.data[
  rownames(coords),"hotspot_module_score1"
]
plot_data <- data.frame(
  tSNE_1=coords[,1],
  tSNE_2=coords[,2],
  hotspot_score=score
)

p <- ggplot(plot_data,
            aes(x=tSNE_1,y=tSNE_2,color=hotspot_score))+
  geom_point(size=1.5,alpha=0.8)+
  scale_color_gradient2(
    low="#3B549D",mid="#F3EED9",high="#B53300",
    midpoint=median(plot_data$hotspot_score)
  )+
  theme_classic()+
  theme(
    axis.line=element_blank(),
    axis.ticks=element_blank(),
    axis.text=element_blank()
  )+
  labs(x=NULL,y=NULL,color="Module score")

print(p)
ggsave("output/Figure4D_scRNA_tsne.svg",
       p,width=12,height=10,units="cm")

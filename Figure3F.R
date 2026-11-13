library(RIdeogram)

####################Figure 3F####################
# 1. 读取已准备的marker数据
markers <- read.csv("data/Figure3F_markers.csv",stringsAsFactors=FALSE)

# 2. 加载RIdeogram内置数据
data(human_karyotype,package="RIdeogram")
data(gene_density,package="RIdeogram")

# 3. 绘图
dir.create("output",showWarnings=FALSE)

ideogram(
  karyotype=human_karyotype,
  overlaid=gene_density,
  label=markers,
  label_type="marker",
  width=170,
  Lx=220,
  Ly=35,
  output="output/Figure3F.svg"
)

# 4. 导出PDF并在RStudio显示
p <- svg2pdf("output/Figure3F.svg",file="output/Figure3F")
grid::grid.newpage()
grid::grid.draw(p)
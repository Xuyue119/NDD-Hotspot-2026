# ============================================================
# Panel A: dRSA 小提琴图
# 输入：
#   <INPUT_DIR>/dRSA_long.csv   —— Gene, Group, Value
#   <INPUT_DIR>/dRSA_stats.csv  —— 单行：Padj_Benign_Non, Padj_Benign_Both,
#                                          Padj_Non_Both, SMD_Non_vs_Both, SMD_CI_L, SMD_CI_U
# 输出：
#   <OUT_DIR>/PanelA_dRSA_violin.pdf
# ============================================================

library(ggplot2)
library(dplyr)
library(ggpubr)

INPUT_DIR <- "<INPUT_DIR>"
OUT_DIR   <- "<OUT_DIR>"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------- 读数据 ----------
drsa_df <- read.csv(file.path(INPUT_DIR, "dRSA_long.csv"),
                    stringsAsFactors = FALSE)
stats   <- read.csv(file.path(INPUT_DIR, "dRSA_stats.csv"),
                    stringsAsFactors = FALSE)

drsa_df$Group <- factor(drsa_df$Group,
                        levels = c("Benign variants",
                                   "Non-clustered variants",
                                   "Both-clustered variants"))

p_bn    <- stats$Padj_Benign_Non[1]
p_bd    <- stats$Padj_Benign_Both[1]
p_nd    <- stats$Padj_Non_Both[1]
smd_nd  <- stats$SMD_Non_vs_Both[1]
smd_l   <- stats$SMD_CI_L[1]
smd_u   <- stats$SMD_CI_U[1]

my_colors <- c(
  "Benign variants"         = "#A1C299",
  "Non-clustered variants"  = "#6C8EBF",
  "Both-clustered variants" = "#D05A4E"
)

# ---------- 显著性比较 ----------
comps <- list()
if (isTRUE(p_bn < 0.05)) comps[[length(comps) + 1]] <- c("Benign variants", "Non-clustered variants")
if (isTRUE(p_bd < 0.05)) comps[[length(comps) + 1]] <- c("Benign variants", "Both-clustered variants")
if (isTRUE(p_nd < 0.05)) comps[[length(comps) + 1]] <- c("Non-clustered variants", "Both-clustered variants")

# ---------- 画图 ----------
p_A <- ggplot(drsa_df, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
  geom_boxplot(width = 0.12, color = "black", alpha = 0.6,
               outlier.shape = NA, linewidth = 0.5) +
  scale_fill_manual(values = my_colors) +
  theme_classic(base_size = 22) +
  labs(
    title    = "A. dRSA (site-local solvent accessibility change magnitude)",
    subtitle = sprintf("SMD(Both\u2212Non) = %+.2f [%+.2f, %+.2f], q = %.3f \u2605",
                       smd_nd, smd_l, smd_u, p_nd),
    x = "", y = "dRSA (|SASA change magnitude|)") +
  theme(
    plot.title    = element_text(size = 24, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 19, color = "#D05A4E", face = "bold", hjust = 0),
    axis.text.x   = element_text(angle = 25, hjust = 1, size = 19, face = "bold"),
    axis.text.y   = element_text(size = 19),
    legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.45)))

if (length(comps) > 0) {
  p_A <- p_A + stat_compare_means(
    comparisons = comps, method = "wilcox.test",
    label = "p.format", hide.ns = TRUE,
    step.increase = 0.14, tip.length = 0.015,
    size = 6.0, vjust = 0.4)
}

ggsave(file.path(OUT_DIR, "PanelA_dRSA_violin.pdf"),
       p_A, width = 12, height = 11)
cat("Panel A saved.\n")
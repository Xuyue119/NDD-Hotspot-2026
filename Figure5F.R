# ============================================================
# Panel B: 7 维度森林图
# 输入：
#   <INPUT_DIR>/panelB_dimension.csv  —— Dimension, SMD, SMD_CI_L, SMD_CI_U, Fisher_Padj
# 输出：
#   <OUT_DIR>/PanelB_Dimension_Forest.pdf
# ============================================================

library(ggplot2)
library(dplyr)

INPUT_DIR <- "<INPUT_DIR>"
OUT_DIR   <- "<OUT_DIR>"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------- 读数据 ----------
df <- read.csv(file.path(INPUT_DIR, "panelB_dimension.csv"),
               stringsAsFactors = FALSE)

df <- df %>%
  mutate(
    Is_Sig     = !is.na(Fisher_Padj) & Fisher_Padj < 0.05,
    Sig_Symbol = ifelse(Is_Sig, "\u2605", "")
  ) %>%
  arrange(SMD)
df$Dimension <- factor(df$Dimension, levels = df$Dimension)

x_range <- range(c(df$SMD_CI_L, df$SMD_CI_U), na.rm = TRUE)
x_span  <- diff(x_range)
x_max   <- x_range[2] + x_span * 0.12
x_min   <- x_range[1] - x_span * 0.05

# ---------- 画图 ----------
p_B <- ggplot(df, aes(x = SMD, y = Dimension, color = Is_Sig)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey30", alpha = 0.6, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = SMD_CI_L, xmax = SMD_CI_U),
                 height = 0.25, linewidth = 1.0) +
  geom_point(aes(shape = Is_Sig, size = Is_Sig)) +
  geom_text(aes(label = Sig_Symbol, x = x_max),
            hjust = 0, size = 10, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(
    values = c("TRUE" = "#D05A4E", "FALSE" = "#888888"),
    breaks = c("TRUE", "FALSE"),
    labels = c("TRUE"  = "\u2605 Fisher q < 0.05",
               "FALSE" = "Fisher q \u2265 0.05"),
    name = "Fisher combination") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
  scale_size_manual(values  = c("TRUE" = 5.0, "FALSE" = 3.5), guide = "none") +
  coord_cartesian(xlim = c(x_min, x_max + x_span * 0.15)) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 19, face = "bold"),
    axis.text.x        = element_text(size = 19),
    plot.title         = element_text(size = 24, face = "bold", hjust = 0),
    plot.subtitle      = element_text(size = 18, color = "grey30", hjust = 0),
    plot.caption       = element_text(size = 15, color = "grey50", hjust = 0),
    legend.position    = "bottom",
    legend.text        = element_text(size = 16),
    legend.title       = element_text(size = 17)) +
  labs(
    x = "SMD (Both \u2212 Non), 95% CI",
    y = "",
    title    = "B. Dimension-level effect sizes and Fisher combination",
    subtitle = "Only solvent exposure reaches Fisher q < 0.05 (\u2605)",
    caption  = "Red filled dot \u2605: Fisher q < 0.05; open dot: q \u2265 0.05. Error bars: 95% CI of SMD.")

ggsave(file.path(OUT_DIR, "PanelB_Dimension_Forest.pdf"),
       p_B, width = 15, height = 9)
cat("Panel B saved.\n")
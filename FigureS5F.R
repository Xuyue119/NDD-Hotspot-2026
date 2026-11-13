# ============================================================
# Panel C: 27 参数森林图
# 输入：
#   <INPUT_DIR>/panelC_parameters.csv  —— Parameter, SMD, SMD_CI_L, SMD_CI_U,
#                                          Padj, Evidence_Level
#                                          Evidence_Level ∈ {Strong, Supporting, No evidence}
# 输出：
#   <OUT_DIR>/PanelC_ParamForest.pdf
# ============================================================

library(ggplot2)
library(dplyr)

INPUT_DIR <- "<INPUT_DIR>"
OUT_DIR   <- "<OUT_DIR>"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------- 读数据 ----------
df <- read.csv(file.path(INPUT_DIR, "panelC_parameters.csv"),
               stringsAsFactors = FALSE)
df <- df %>%
  filter(!is.na(SMD)) %>%
  mutate(
    Evidence_Level = factor(Evidence_Level,
                            levels = c("Strong", "Supporting", "No evidence")),
    Param_f = factor(Parameter, levels = Parameter[order(SMD)])
  )

# ---------- 画图 ----------
p_C <- ggplot(df, aes(x = SMD, y = Param_f, color = Evidence_Level)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_errorbarh(aes(xmin = SMD_CI_L, xmax = SMD_CI_U),
                 height = 0.3, linewidth = 0.7) +
  geom_point(size = 4.0) +
  scale_color_manual(
    values = c("Strong" = "#D05A4E", "Supporting" = "#E8A33D", "No evidence" = "#BDBDBD"),
    labels = c("Strong" = "FDR-significant (q < 0.05)",
               "Supporting" = "CI excl. 0, q \u2265 0.05",
               "No evidence" = "No evidence")) +
  theme_bw(base_size = 20) +
  labs(x = "SMD (Both \u2212 Non)", y = "",
       title = "C. All 27 non-redundant parameters",
       subtitle = "dRSA is the only FDR-significant parameter (red)",
       color = "Evidence level") +
  theme(axis.text.y   = element_text(size = 15),
        axis.text.x   = element_text(size = 17),
        plot.title    = element_text(size = 24, face = "bold"),
        plot.subtitle = element_text(size = 18, color = "grey40"),
        legend.position = "bottom",
        legend.text     = element_text(size = 16),
        legend.title    = element_text(size = 17))

ggsave(file.path(OUT_DIR, "PanelC_ParamForest.pdf"),
       p_C, width = 15, height = 16)
cat("Panel C saved.\n")
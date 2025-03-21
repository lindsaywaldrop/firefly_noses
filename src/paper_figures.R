# Creating paneled figures for paper
library(cowplot)
library(ggplot2)
library(patchwork)

# Figure 1

p_sem1 <- ggdraw() + draw_image("./doc/figures/figure1/p_macdermotti_antennae.tif")
p_sem2 <- ggdraw() + draw_image("./doc/figures/figure1/Lucidota_sp_192_f1W.tif")
p_sems_fig <- p_sem1 + p_sem2 + plot_layout(widths = c(1, 1.63)) +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(hjust = -0.5, vjust = 2.5, size = 20,
                                color = "white", face = "plain"))
ggsave("./doc/figures/figure1/SEM_fig.png", p_sems_fig, height = 5, width = 10, dpi = 600)


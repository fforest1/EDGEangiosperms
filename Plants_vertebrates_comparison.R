# Comparison of threatened evolutionary history, EDGE species numbers and threatened species
# numbers across selected groups of plants and animals

library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggimage)
library(rsvg)
library(grid)
library(tidyr)

########## Panel A - boxplots threatened evolutionary history for each clade

perc_ePDloss <- read.csv(file="perc_thrPD.csv")

plant_clades <- c("Angiosperms", "Gymnosperms")

ordered_clades <- perc_ePDloss %>%
  filter(!clade %in% plant_clades) %>%
  group_by(clade) %>%
  summarize(mean_loss = mean(perc.ePDloss, na.rm = TRUE)) %>%
  arrange(desc(mean_loss)) %>%
  pull(clade)

order <- c(plant_clades, ordered_clades)

perc_ePDloss$clade <- factor(perc_ePDloss$clade, levels = order)

perc_ePDloss$color_group <- ifelse(perc_ePDloss$clade %in% plant_clades, "plants", "vertebrates")

pA <- 
ggplot(perc_ePDloss, aes(x = clade, y = perc.ePDloss, fill = color_group, color = color_group)) +
  geom_boxplot(alpha = 0.5, size = 0.1, outlier.size = 1) +
  scale_fill_manual(values = c("plants" = "darkgreen", "vertebrates" = "blue")) +
  scale_color_manual(values = c("plants" = "darkgreen", "vertebrates" = "blue")) +
  labs(x = "Clades", y = "% Threatened evolutionray history") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "grey"), 
        axis.ticks.y = element_line(colour = "grey"))

########## Panel B - Proportion (in %) of EDGE species and threatened species in all major groups of plants and animals 
clades <- read.csv(file = "clade comparison.csv")
clades$Clade <- factor(clades$Clade, levels = clades$Clade) 

clades_long <- clades %>%
  pivot_longer(cols = c(perc.edge.spp, perc.thr),
               names_to = "Metric",
               values_to = "Value")
pB <-
  ggplot(clades_long, aes(x = Clade, y = Value, fill = Clade)) +
  geom_col(aes(group = Metric, alpha = Metric), position = position_dodge(width = 0.8), width = 0.8, color = "black", linewidth = 0.1) +
  scale_fill_manual(values = c(rep("darkgreen", 2), rep("dodgerblue3", 9)), guide = "none") +
  scale_alpha_manual(values = c("perc.edge.spp" = 1, "perc.thr" = 0.5)) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 0.5) + 
  labs(x = "Categories", y = "% Species") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(colour = "grey"),  
    axis.ticks.y = element_line(colour = "grey"),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, max(clades_long$Value) * 1.2))

pA <- pA + theme(plot.margin = margin(t = 10, b = 10, unit = "pt"))
pB <- pB + theme(plot.margin = margin(t = 10, b = 10, unit = "pt"))

base_text_size <- 7  

pA <- pA + theme(
  text = element_text(size = base_text_size),         
  axis.text = element_text(size = base_text_size),
  axis.title = element_text(size = base_text_size),
  legend.text = element_text(size = base_text_size),
  legend.title = element_text(size = base_text_size),
  strip.text = element_text(size = base_text_size)
)

pB <- pB + theme(
  text = element_text(size = base_text_size),
  axis.text = element_text(size = base_text_size),
  axis.title = element_text(size = base_text_size),
  legend.text = element_text(size = base_text_size),
  legend.title = element_text(size = base_text_size),
  strip.text = element_text(size = base_text_size)
)

g <- arrangeGrob(pA, pB, layout_matrix = matrix(c(1, 2, 3), ncol = 1))
ggsave("Figure1_new.pdf", plot = g, units = "cm", width = 9, height = 17, dpi = 300)

# With minor manual cosmetic alterations. Icons were added manually. 

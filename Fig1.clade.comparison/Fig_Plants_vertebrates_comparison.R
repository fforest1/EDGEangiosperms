# Figure 1 - Comparison of threatened evolutionary history and EDGE species numbers across selected groups of plants and animals

library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggimage)
library(rsvg)
library(grid)

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
  geom_boxplot(alpha = 0.5, size = 0.2, outlier.size = 1) +
  scale_fill_manual(values = c("plants" = "darkgreen", "vertebrates" = "blue")) +
  scale_color_manual(values = c("plants" = "darkgreen", "vertebrates" = "blue")) +
  labs(x = "Clades", y = "% Threatened evolutionray history") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "grey"), 
        axis.ticks.y = element_line(colour = "grey"))


########## Panel B - threatened evolutionary history per species for each clade against their total species richness  (log)
clades <- read.csv(file = "clade comparison.csv")

pB <- 
ggplot(clades, aes(x = total.spp, y = evo.loss.per.spp, color = group)) +
  geom_point(size = 2) +
  geom_text(aes(label = Clade), vjust = -0.5, size = 3) +  # Adding labels
  labs(x = "Species richness (log)", y = "Threatened evolutionary history per species") +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000, 100000, 340000),
    labels = scales::comma)  +
  scale_color_manual(values = c("plants" = "darkgreen", "vertebrates" = "dodgerblue3")) +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "grey"), 
        axis.ticks.y = element_line(colour = "grey"))


########## Panel C - Proportion (in %) of EDGE species in all major groups of plants and animals 
#clades$Clade <- factor(clades$Clade, levels = unique(clades$Clade))
clades$Clade <- factor(clades$Clade, levels = clades$Clade) 

pC <- 
  ggplot(clades, aes(x = Clade, y = perc.edge.spp, fill = Clade)) +
  geom_col(color = "black",  # Border color of the bars
           linewidth = 0.3,  # Adjusts the border width
           alpha = 0.5) +  # Set transparency for the bars (adjust as needed)
  scale_fill_manual(values = c(rep("darkgreen", 2), rep("dodgerblue3", 9))) +  # Set fill color for each group
  labs(x = "Categories", y = "% EDGE species") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Adjust x-axis labels angle and size
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),  # Remove grid lines
    legend.position = "none",
    axis.line.x = element_blank(),  # Removes x-axis line
    axis.line.y = element_line(colour = "grey"),  
    axis.ticks.y = element_line(colour = "grey")
    )+
  coord_cartesian(ylim = c(0, max(clades$perc.edge.spp) * 1.2))


# Add extra margin space to each plot
pA <- pA + theme(plot.margin = margin(t = 10, b = 10, unit = "pt"))
pB <- pB + theme(plot.margin = margin(t = 10, b = 10, unit = "pt"))
pC <- pC + theme(plot.margin = margin(t = 10, b = 10, unit = "pt"))

# change text size throughout
base_text_size <- 7  

pA <- pA + theme(
  text = element_text(size = base_text_size),         # Global text size
  axis.text = element_text(size = base_text_size),    # Axis tick labels
  axis.title = element_text(size = base_text_size),   # Axis titles
  legend.text = element_text(size = base_text_size),  # Legend text
  legend.title = element_text(size = base_text_size), # Legend title
  strip.text = element_text(size = base_text_size)    # Facet labels
)

pB <- pB + theme(
  text = element_text(size = base_text_size),
  axis.text = element_text(size = base_text_size),
  axis.title = element_text(size = base_text_size),
  legend.text = element_text(size = base_text_size),
  legend.title = element_text(size = base_text_size),
  strip.text = element_text(size = base_text_size)
)

pC <- pC + theme(
  text = element_text(size = base_text_size),
  axis.text = element_text(size = base_text_size),
  axis.title = element_text(size = base_text_size),
  legend.text = element_text(size = base_text_size),
  legend.title = element_text(size = base_text_size),
  strip.text = element_text(size = base_text_size)
)

g <- arrangeGrob(pA, pB, pC, layout_matrix = matrix(c(1, 2, 3), ncol = 1))

# Save the figure with spacing
ggsave("Figure1e.pdf", plot = g, units = "cm", width = 9, height = 17, dpi = 300)

# Icons were added manually. 

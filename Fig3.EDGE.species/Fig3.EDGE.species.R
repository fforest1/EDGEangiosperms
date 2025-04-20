# Figure 3 - EDGE species

library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)

options(scipen = 999)
edge.scores <- read.csv(file="Forest_etal_SupplMat_S10_200EDGEscores.csv")
edge.scores <- edge.scores[, -1]
edge.scores.t <- pivot_longer(edge.scores, cols = -Species.1, values_to = "EDGE.score")

top_25 <- edge.scores.t %>%
  group_by(Species.1) %>%
  summarize(median_value = median(EDGE.score, na.rm = TRUE)) %>%
  arrange(desc(median_value)) %>%
  slice_head(n = 25)

top_25 <- edge.scores.t %>% filter(Species.1 %in% top_25$Species.1)

edge.summary <- read.csv(file="Forest_etal_SupplMat_S1_EDGEspp_RL&Pred.csv")
edge.summary <- edge.summary[, c("Species", "threat")]
top_25 <- merge(top_25, edge.summary, by.x="Species.1", by.y="Species", all.x=TRUE)

top_25$Species.1 <- gsub("_", " ", top_25$Species.1)

pA <- 
ggplot(top_25, aes(x = reorder(Species.1, EDGE.score, median), y = EDGE.score, fill = threat)) +
  geom_boxplot(outlier.shape = NA, size = 0.1) +
  labs(x="", y = "EDGE score (MY at risk)") +
  ylim(0, 100) +
  scale_fill_manual(values = c("CR" = "red", "EN" = "orange", "thr" = "darkorchid")) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.y = element_text(face = "italic", size = 6),
        axis.line.x = element_line(color = "black", size = 0.1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(color = "black", size = 0.1),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        panel.grid.major.x = element_line(color = "grey80", linewidth = 0.2, linetype = "dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()) +
   guides(fill = guide_legend(title = NULL)) +
  coord_flip() 

ggsave("Figure3.pdf", width = 18, height = 12, units = "cm", limitsize = FALSE)

# minor adjustments and addition of images done manually

# Figure 5 - EDGE families

library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)

options(scipen = 999)
edge.summary <- read.csv(file="Forest_etal_SupplMat_S1_EDGEspp_RL&Pred.csv")
edge.summary <- edge.summary[, c("Species", "Family")]

edge.scores <- read.csv(file="Forest_etal_SupplMat_S10_200EDGEscores.csv")
edge.scores <- edge.scores[, -1]
edge.scores.t <- pivot_longer(edge.scores, cols = -Species.1, values_to = "EDGE.score")

edge.scores.t <- merge(edge.scores.t, edge.summary, by.x="Species.1", by.y="Species", all.x=TRUE)

edge.fam <- read.csv(file="Forest_etal_SupplMat_S3_EDGEfam_RL&Pred.csv")
edge.fam <- edge.fam[edge.fam$is.EDGE.fam == 'Yes',] # only EDGE families
edge.fam <- edge.fam[order(edge.fam$N.edgespp, decreasing = TRUE), ] # rank by N.edgespp
edge.fam <- edge.fam[1:25, ] # keep top 25

edge.scores.t <- edge.scores.t %>% semi_join(edge.fam, by = "Family")  # trim edge.scores.t to only species from top 25 EDGE families
edge.scores.t$name <- NULL 

# Aggregate mean EDGE.score by Species.1
edge.scores.t.mean <- edge.scores.t %>%
  group_by(Species.1, Family) %>%
  summarise(mean_EDGE = mean(EDGE.score, na.rm = TRUE), .groups = "drop")

edge.fam <- edge.fam[, c("Family", "N.edgespp", "P.thrt")]
merged_data <- merge(edge.scores.t.mean, edge.fam, by = "Family")

# Reorder 'family' based on 'N.edgespp'
merged_data$Family <- reorder(merged_data$Family, merged_data$N.edgespp, FUN = median, decreasing = FALSE)

# Create the boxplot with the reordered 'family'
pB <- 
ggplot(merged_data, aes(x = Family, y = mean_EDGE, fill = P.thrt)) +
  geom_boxplot(outlier.shape = NA, size = 0.1) +
  labs(x="", y = "Mean EDGE score (MY at risk)") +
  ylim(1, 32) +
  scale_fill_gradient(low = "purple", high = "yellow") +  
  theme_minimal() +
  theme(legend.position = "right",
        axis.line.x = element_line(color = "black", size = 0.1),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.x = element_line(color = "black", size = 0.1),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        panel.grid.major.x = element_line(color = "grey80", linewidth = 0.2, linetype = "dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  coord_flip()

ggsave("Figure5.pdf", width = 18, height = 12, units = "cm", limitsize = FALSE)

# minor adjustments and addition of images done manually

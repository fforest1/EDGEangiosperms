# Figure 2 - Distribution of EDGE species and EDGE families across the tree of life of angiosperms

library(ggplot2)
library(ape)
library(ggtree)
library(treeio)
library(ggtreeExtra)
library(dplyr)

tr_fam <- read.tree(file = "family_tree.tre")

edge_fam <- read.csv(file="Forest_etal_SupplMat_S3_EDGEfam_RL&Pred.csv")
med.edge <- edge_fam[, c("Family", "edge.med")]

# restrict EDGE families to 200 spp
edge_fam$N.edgespp2 <- ifelse(edge_fam$N.edgespp > 200, 200, edge_fam$N.edgespp)

# Function to get MRCA or handle single-tip groups
get_mrca_safe <- function(tr_fam2, Family) {
  if (length(Family) == 1) {
    return(which(tr_fam2$tip.label == Family))  # Return the tip index if only one tip
  } else {
    return(getMRCA(tr_fam2, Family))  # Use getMRCA() for multiple tips
  }
}

# Find MRCA (or tip index for single-tip groups)
mrca_nodes <- edge_fam %>% group_by(order) %>% summarize(node = get_mrca_safe(tr_fam, Family), .groups = "drop")  
mrca_nodes$node <- as.numeric(mrca_nodes$node)

p <- 
  ggtree(tr_fam, ladderize = TRUE, size=0.2, layout="fan", open.angle=10, aes(color = pmin(pmax(edge.med, 0.1), 38))) %<+% med.edge +
  
  scale_color_gradientn(colors = c("darkblue", "blue", "cyan", "yellow", "orange", "orange", "red", "darkred"), 
                        trans = "log",  
                        breaks = c(0.1, 0.2, 1, 5, 38),
                        limits = c(0.05, 40),
                        labels = scales::label_number(),
                        name = "Median family EDGE score (myr)") +  

  geom_fruit(data=edge_fam, geom=geom_bar,
             mapping=aes(y=Family, x=N.edgespp2, fill=is.EDGE.fam),
             pwidth=0.38, 
             orientation="y", 
             stat="identity",
             offset = 0.35,
             axis.params=list(
               axis       = "x",
               text.size  = 1.5,
               hjust      = 0.5,
               vjust      = 1.5,
               nbreak     = 3,
             ),
             grid.params=list()
  ) + 
  
  scale_fill_manual(values=c("#696969","#FF0000"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.3, order=4),
                    name = "EDGE families") +
  
  annotate("text", x = 210, y = 0.1, label = "No EDGE species", size = 2) +
  
  theme(legend.position.inside = c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  )

# add order labels
p + geom_cladelab(data = mrca_nodes, 
                       mapping = aes(node = node, label = order), 
                       align = TRUE,  
                       offset = 1.5,  
                       offset.text = 1,
                       fontsize = 1, 
                       barsize = 1,
                       barcolour = "#696969",
                       angle = "auto")  


ggsave("Figure2.pdf", width = 18, height = 26, units = "cm", limitsize = FALSE)

# minor adjustments made separately.

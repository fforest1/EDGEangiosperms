######### Modified EDGE2 calculation code - for adding a category-based nudge to pext scores ######### 

# Original code by Rikki Gumbs
# Modified by Ruth Brown rcb22@ic.ac.uk
# Adapted by Matilda Brown m.brown2@kew.org

# Provide phylogenetic tree and dataframe with Red List categories and predictions for 200 draws: 
# function returns three objects: 
# 1. dataframe with terminal branch length, GE2, ED2 and EDGE2 scores for each species
# 2. expected PD loss tree
# 3. PD and expected PD loss in Myr for the clade

library(ape)
library(dplyr)

tree <- read.tree(file="Forest_etal_SupplMat_S8_200trees.tre")

#To avoid stochastic differences that may arise from inter-draw differences, we simulated an increase in GE2 score directly, without resampling from the curve. 
pext.dist <- read.csv("pext_dist.csv")

# Adjustement value; we used 0.2, 0.5 and 1
adjustment_amount <- 0.2

## FUNCTION DEFINITIONS ----
# score adjust function
score_adjust <- function(initial_scores, calibration_data = pext.dist, category_increase){

initial_scores <- initial_scores[order(initial_scores)]

initial_cat_scores <-  approx(x = calibration_data$pext, 
                              y = calibration_data$RL.cat.num, 
                              xout = initial_scores, ties="ordered")$y

 df <- data.frame(initialGE=initial_scores,
               initialRL=initial_cat_scores) %>% 
      mutate(newRL = initialRL + category_increase) %>% 
      mutate(newRL = case_when(
        newRL>5~5,
        newRL<0~0,
        TRUE~newRL
      ))
    
    newGEscores <- approx(x = calibration_data$RL.cat.num, 
                            y = calibration_data$pext, 
                            xout = df$newRL, ties="ordered")$y
    df %>% 
      mutate(newGE = newGEscores) %>% 
      mutate(newGE = case_when(
        initialGE == 1 ~ 1,
        newGE>0.9999 ~ 0.9999, 
        newGE<0.0001 ~ 0.0001,
        TRUE ~ newGE
      ))
}

# remove excess pext, and reorders to same order as tree$tip.label
into_order <- function(tree, pext){
  new_pext <- pext[match(tree$tip.label, pext$species),]
return (new_pext)
}

# order tree components
reorder_tree <- function(tree, ordering){
  tree@edge.length <- tree@edge.length[ordering]
  tree@edge <- tree@edge[ordering,]
  return(tree)
}

# EDGE2 Function
EDGE2_mod <- function(tree, pext){
  
  require(phylobase)
  require(data.table)
  
  names(pext) <- c("species","pext")
  
  N_species <- length(tree$tip.label)
  N_nodes <- tree$Nnode
  N_tot <- N_species + N_nodes
  
  # ensure extinction probabilities are given in same order as tree$tip.label.
  if (!identical(tree$tip.label, pext$species)){
    pext <- into_order(tree, pext)
  }
  
  if(!class(tree) == "phylo"){
    tree <- as(tree, "phylo")
  }
  
  tree_dat <- data.frame(Species = as.character(tree$tip.label),
                         TBL = NA, 
                         pext = pext$pext, ED = NA, EDGE = NA)
  ePD.dat <- data.frame(PD = sum(tree$edge.length),ePDloss = NA)
  
  tree <- as(tree, "phylo4")
  root <- rootNode(tree)
  nodes <- c(root, descendants(tree, root, "all"))
  
  # reorder tree components more instinctively, such that nodes are easier to find
  ord <- order(nodes)
  tree <- reorder_tree(tree, ord)
  nodes <- nodes[ord]
  
  tree_dat$TBL <- tree@edge.length[1:N_species]
  
  node_data <- data.frame(Node = 1:N_tot, Pext = rep(1, N_tot), Edge_Sum = NA)
  node_data[1:length(tree@label), 2] <- pext[,2]
  
  # assign the product of its descendant tips to each node
  for (i in c(1:length(tree@label), N_tot:(root+1))){         # for each node, beginning with tips
      anc <- tree@edge[i,1]                                   # find ancestor of node
      node_data[anc, 2] <- node_data[anc, 2]*node_data[i,2]   # muliply ancestor value by node "pext" 
  }

  # multiply each edge.length by each pext caluclated above
  for(i in 1:length(nodes)){
    tree@edge.length[i] <- tree@edge.length[i]*node_data[i,2]
  }
  save(tree, file ="tree.rda")

  if (is.na(tree@edge.length[root])){
    tree@edge.length[root] <- 0
  }
  node_data$Edge_Sum[root] <- tree@edge.length[root]
  
  # for each internal node, summate ancesteral edgelengths
  for (i in (root+1):N_tot){
    ans <- tree@edge[i,1]
    node_data$Edge_Sum[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
  }
  
  # for each tip, summate ancesteral edgelengths to find EDGE2 score
  for (i in 1:N_species){
    ans <- tree@edge[i,1]
    tree_dat$EDGE[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
  }  
  
  tree_dat$ED <- tree_dat$EDGE / tree_dat$pext
  # reorder tree
  tree <- reorder_tree(tree, order(ord))
 
  tree <- as(tree, "phylo")
  ePD.dat$ePDloss <- sum(tree$edge.length)
  edge.res <- list(tree_dat,tree,ePD.dat)
  return(edge.res)
}

######### END FUNCTIONS ---- 

######### EDGE2 CALC (in loop) ----  
EDGE.2.list <- list()
exp.PD.list <- list()
exp.PD.trees <- list()

for(i in 1:length(tree)){

Species <- read.csv(file = paste("pext_dist_1/pext_tree",i,".csv", sep = ""))

Species$id <- 1:nrow(Species) 
Species <- Species[order(Species$pext),]
Species$pext <- score_adjust(Species$pext, pext.dist, adjustment_amount)$newGE
Species <- Species[order(Species$id),]
Species <- Species[,c("Species","RL.cat","pext")]

print((paste(i, "-1")))

  #rename and order columns for function
  col.num.1 <- which(colnames(Species) == "Species")
  col.num.2 <- which(colnames(Species) == "pext")
  sp.pext <- Species[,c(col.num.1,col.num.2)]
  names(sp.pext) <- c("Species","pext")
print((paste(i, "-2")))

  # calculate EDGE 2
  res  <- EDGE2_mod(tree[[i]],sp.pext)
print((paste(i, "-3")))

  # determine whether species are above median EDGE for each iteration for ease later
  res2 <- data.frame(res[[1]],above.median = 0, Iteration = i)
  res2$above.median[res2$EDGE > median(res2$EDGE)] <- 1

  # put ED/EDGE scores in list 
  EDGE.2.list[[length(EDGE.2.list)+1]] <- res2[order(res2$EDGE,decreasing = T),]

  #calculate PD, ePDloss and ePD for the iteration and put in list
  exp.PD.list[[length(exp.PD.list)+1]] <- data.frame(PD = sum(tree[[i]]$edge.length),ePDloss = sum(res[[2]]$edge.length),
                                                  ePD = (sum(tree[[i]]$edge.length) - sum(res[[2]]$edge.length)))
  
  # put the pext-transformed ePDloss tree in list
  exp.PD.trees[[length(exp.PD.trees)+1]] <- res[[2]]
print((paste(i, "-end")))
}

## WRITE OUTPUT----
save.image(file="EDGE2_increasing_pext02.results")



######### Modified EDGE2 calculation code ######### 

# Original code by Rikki Gumbs
# Modified by Ruth Brown
# Adapted to run with predictions by Rikki Gumbs and Felix Forest

# Provide phylogenetic tree and dataframe with Red List categories and predictions for 200 draws: 
# function returns three objects: 
# 1. dataframe with terminal branch length, GE2, ED2 and EDGE2 scores for each species
# 2. expected PD loss tree
# 3. PD and expected PD loss in Myr for the clade

library(ape)

tree <- read.tree(file="Forest_etal_SupplMat_S8_200trees.tre")
pred <- read.csv(file="Forest_etal_SupplMat_S9_Pred_200draws.csv")

pext <- rev(c(0.97, 0.97/2, 0.97/4,0.97/8,0.97/16))
GE.2.calc <- function(pext){
  require(geiger)
  treesim <- sim.bdtree(n=10000)
  iucn <- sample(1:5, size=length(treesim$tip.label), replace=TRUE)
  data <- data.frame(species=treesim$tip.label, pext=pext[iucn])
  data <- data[order(data$pext),]
  data$rank <- seq_len(nrow(data))
  rank <- c(0, with(data, tapply(rank, pext, median)))
  pext <- c(0, pext)
  rank.sq <- rank^2; rank.cub <- rank^3; rank.qu <- rank^4; rank.quu <- rank^5
  model <- lm(pext ~ rank + rank.sq + rank.cub + rank.qu)
  data$rank.sq <- data$rank^2; data$rank.cub <- data$rank^3; data$rank.qu <- data$rank^4; data$rank.quu <- data$rank^5
  data$rank.pext <- predict(model, data)
  data$rank.pext[data$rank.pext <= 0] <- 0.0001
  data$rank.pext[data$rank.pext >= 1] <- 0.9999
  pext.LC <- data.frame(RL.cat = "LC", pext =data$rank.pext[data$pext == pext[2]])
  pext.NT <- data.frame(RL.cat = "NT", pext =data$rank.pext[data$pext == pext[3]])
  pext.VU <- data.frame(RL.cat = "VU", pext =data$rank.pext[data$pext == pext[4]])
  pext.EN <- data.frame(RL.cat = "EN", pext =data$rank.pext[data$pext == pext[5]])
  pext.CR <- data.frame(RL.cat = "CR", pext =data$rank.pext[data$pext == pext[6]])
  return(rbind(pext.CR,pext.EN, pext.VU, pext.NT, pext.LC))
}
pext.dist <- GE.2.calc(pext)

EDGE.2.list <- list()
exp.PD.list <- list()
exp.PD.trees <- list()

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


for(i in 1:length(tree)){
  j <- i + 1
  Species <- pred[, c(1, j)]
  names(Species)[names(Species) == paste("draw", i, sep="")] <- "RL.cat"
  Species$RL.cat <- gsub("0","not", Species$RL.cat)
  Species$RL.cat <- gsub("1","thr", Species$RL.cat)

  # To run with Red List categories only (omitting predictions), add "#" on lines 144-145 and remove "#" from lines 148-149.
  # Species$RL.cat <- gsub("0","NE", Species$RL.cat)
  # Species$RL.cat <- gsub("1","NE", Species$RL.cat)
  
  # randomly select pext scores for all species in tree based on their GE
  Species$pext <- 0
  Species$pext[which(Species$RL.cat == "LC")] <- sample(pext.dist$pext[pext.dist$RL.cat == "LC"],
                                                        length(Species$RL.cat[which(Species$RL.cat == "LC")]),replace = T)
  Species$pext[which(Species$RL.cat == "NT")] <- sample(pext.dist$pext[pext.dist$RL.cat == "NT"],
                                                        length(Species$RL.cat[which(Species$RL.cat == "NT")]),replace = T)
  Species$pext[which(Species$RL.cat == "VU")] <- sample(pext.dist$pext[pext.dist$RL.cat == "VU"],
                                                        length(Species$RL.cat[which(Species$RL.cat == "VU")]),replace = T)
  Species$pext[which(Species$RL.cat == "EN")] <- sample(pext.dist$pext[pext.dist$RL.cat == "EN"],
                                                        length(Species$RL.cat[which(Species$RL.cat == "EN")]),replace = T)
  Species$pext[which(Species$RL.cat %in% c("CR","EW"))] <- sample(pext.dist$pext[pext.dist$RL.cat == "CR"],
                                                                 length(Species$RL.cat[which(Species$RL.cat == "CR")]),replace = T)
  Species$pext[which(Species$RL.cat == "thr")] <- sample(pext.dist$pext[pext.dist$RL.cat %in% c("VU", "EN", "CR")],
                                                         length(Species$RL.cat[which(Species$RL.cat == "thr")]),replace = T)
  Species$pext[which(Species$RL.cat == "not")] <- sample(pext.dist$pext[pext.dist$RL.cat %in% c("LC", "NT")],
                                                         length(Species$RL.cat[which(Species$RL.cat == "not")]),replace = T)
  Species$pext[which(Species$RL.cat %in% c("EX"))] <- 1

  # To run with Red List categories only (omitting predictions), add "#" on lines 163-166 and remove "#" from lines 171-172.
  # NE (and DD) selected from entire GE2 distribution 
  # Species$pext[which(Species$RL.cat == "NE")] <- sample(pext.dist$pext[pext.dist$RL.cat %in% c("VU", "EN", "CR", "LC", "NT")],
  #                                                       length(Species$RL.cat[which(Species$RL.cat == "NE")]),replace = T)
  
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


save.image(file="EDGE2.results")


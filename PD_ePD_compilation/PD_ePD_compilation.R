
# Compilation of total and threatened evolutionary history

library(ape)

##################################################
# compile total evolutionary history for 200 trees
trees <- read.tree(file="Forest_etal_EDGEangio_DataS4_200trees.tre")

total.PD <- numeric(length(trees))  

for(i in seq_along(trees)) {
  total.PD[i] <- sum(trees[[i]]$edge.length)
  print(i)
}

median(total.PD)


##################################################
# compile threatened evolutionary history for 200 weighted trees
w.trees <- read.tree(file="Forest_etal_DataS6_200_Weighted_trees.tre")

thr.ePD <- numeric(length(w.trees))  

for(i in seq_along(w.trees)) {
  thr.ePD[i] <- sum(w.trees[[i]]$edge.length)
  print(i)
}

median(thr.ePD)




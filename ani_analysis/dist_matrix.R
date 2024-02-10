#ANI analysis
library(dplyr)
library(tidyr)
library(ape)
#mash distances matrix
distances <- read.delim("F:/zyf/mbdata/prok/distances.tab", header=FALSE)
distances <- distances[-5]
distance_matrix_df <- distances %>%
  select(V1, V2, V4) %>% 
  pivot_wider(names_from = V2, values_from = V4, values_fill = list(V4 = 0)) 
# create a new column for each unique value in V2, taking values from V4

distance_matrix <- as.matrix(distance_matrix_df)
distance_matrix <- distance_matrix[,-1]
rownames(distance_matrix) <- colnames(distance_matrix)

#create and plot NJ tree
distance_vector <- as.dist(distance_matrix)
nj_tree <- bionj(distance_vector)
plot(nj_tree, show.tip.label = FALSE)
write.csv(distance_matrix, file = "distance_matrix.csv", row.names = TRUE, col.names = TRUE)
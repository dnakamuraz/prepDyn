#!/usr/bin/env Rscript

#############
# PRE-AMBLE #
#############

# Set the working directory 
setwd("/Users/labanfibios/Desktop/Doutorado/Project/B3_PrepDyn/GitHub/test_data/empirical_data/")

# Load required packages
library(gplots)
library(pheatmap)
if (!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)
library(RNODE)
library(ggplotify)
library(patchwork)

#############
# LOAD DATA #
#############

# Empirical 1 (de Sá)
# Load consensus
e1_1_noPrep = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/1_noPreprocessing/v1.4_TNT_strictConsensus.nwk")
e1_2_noPart = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/2_noPartitioning/v2.4_TNT_strictConsensus.nwk")
e1_3_max = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/3_max/v6.4_TNT_strictConsensus.nwk")
e1_4_balanced = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/4_balanced/round1/v3.1.4_TNT_strictConsensus.nwk")
e1_5_conservative = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/5_conservative/round1/v4.1.4_TNT_strictConsensus.nwk")
e1_6_equal = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/6_equal/round1/v5.1.4_TNT_strictConsensus.nwk")
# Create a list of consensus trees 1
e1 = c(e1_1_noPrep, e1_2_noPart, e1_3_max, e1_4_balanced, e1_5_conservative, e1_6_equal)
# Load MPTs
mpt1_1_noPrep = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/1_noPreprocessing/v1.3_PHYG_Sa_bestUnique.nwk")
mpt1_2_noPart = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/2_noPartitioning/v2.3_PHYG_Sa_bestUnique.nwk")
mpt1_3_max = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/3_max/v6.3_PHYG_Sa_bestUnique.nwk")
mpt1_4_balanced = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/4_balanced/round1/v3.1.3_PHYG_Sa_bestUnique.nwk")
mpt1_5_conservative = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/5_conservative/round1/v4.1.3_PHYG_Sa_bestUnique.nwk")
mpt1_6_equal = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical1_deSa/6_equal/round1/v5.1.3_PHYG_Sa_bestUnique.nwk")
# Create a list of MPTs 1
mpt1 = c(mpt1_1_noPrep, mpt1_2_noPart, mpt1_3_max, mpt1_4_balanced, mpt1_5_conservative, mpt1_6_equal)

# Empirical 2 (Chiari)
e2_1_noPrep = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/1_noPreprocessing/v1.4_TNT_strictConsensus.nwk")
e2_2_noPart = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/2_noPartitioning/v2.4_TNT_strictConsensus.nwk")
e2_3_max = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/3_max/v6.4_TNT_strictConsensus.nwk")
e2_4_balanced = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/4_balanced/round4/v3.4.4_TNT_strictConsensus.nwk")
e2_5_conservative = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/5_conservative/round2/v4.2.4_TNT_strictConsensus.nwk")
e2_6_equal = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/6_equal/round4/v5.4.4_TNT_strictConsensus.nwk")
# Create a list of consensus trees 2
e2 = c(e2_1_noPrep, e2_2_noPart, e2_3_max, e2_4_balanced, e2_5_conservative, e2_6_equal)
# Load MPTs
mpt2_1_noPrep = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/1_noPreprocessing/v1.3_PHYG_Chiari_bestUnique.nwk")
mpt2_2_noPart = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/2_noPartitioning/v2.3_PHYG_Chiari_bestUnique.nwk")
mpt2_3_max = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/3_max/v6.3_PHYG_Chiari_bestUnique.nwk")
mpt2_4_balanced = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/4_balanced/round4/v3.4.3_PHYG_Chiari_bestUnique.nwk")
mpt2_5_conservative = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/5_conservative/round2/v4.2.3_PHYG_Chiari_bestUnique.nwk")
mpt2_6_equal = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical2_Chiari/6_equal/round4/v5.4.3_PHYG_Chiari_bestUnique.nwk")
# Create a list of MPTs 2
mpt2 = c(mpt2_1_noPrep, mpt2_2_noPart, mpt2_3_max, mpt2_4_balanced, mpt2_5_conservative, mpt2_6_equal)

# Empirical 3 (Nakamura)
e3_1_noPrep = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/1_noPreprocessing/v1.4_TNT_strictConsensus.nwk")
e3_2_noPart = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/2_noPartitioning/v2.4_TNT_strictConsensus.nwk")
e3_3_max = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/3_max/v6.4_TNT_strictConsensus.nwk")
e3_4_balanced = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/4_balanced/round3/v3.3.4_TNT_strictConsensus.nwk")
e3_5_conservative = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/5_conservative/round2/v4.2.4_TNT_strictConsensus.nwk")
e3_6_equal = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/6_equal/round5/v5.5.4_TNT_strictConsensus.nwk")
# Create a list of consensus trees 3
e3 = c(e3_1_noPrep, e3_2_noPart, e3_3_max, e3_4_balanced, e3_5_conservative, e3_6_equal)
# Load MPTs
mpt3_1_noPrep = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/1_noPreprocessing/v1.3_PHYG_Nakamura_bestUnique.nwk")
mpt3_2_noPart = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/2_noPartitioning/v2.3_PHYG_Nakamura_bestUnique.nwk")
mpt3_3_max = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/3_max/v6.3_PHYG_Nakamura_bestUnique.nwk")
mpt3_4_balanced = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/4_balanced/round3/v3.3.3_PHYG_Nakamura_bestUnique.nwk")
mpt3_5_conservative = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/5_conservative/round2/v4.2.3_PHYG_Nakamura_bestUnique.nwk")
mpt3_6_equal = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical3_Nakamura2/6_equal/round5/v5.5.3_PHYG_Nakamura_bestUnique.nwk")
# Create a list of MPTs 3
mpt3 = c(mpt3_1_noPrep, mpt3_2_noPart, mpt3_3_max, mpt3_4_balanced, mpt3_5_conservative, mpt3_6_equal)

# Empirical 4 (Domburg)
e4_1_noPrep = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/1_noPreprocessing/v1.4_TNT_strictConsensus.nwk")
e4_2_noPart = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/2_noPartitioning/v2.4_TNT_strictConsensus.nwk")
e4_3_max = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/3_max/v6.4_TNT_strictConsensus.nwk")
e4_4_balanced = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/4_balanced/round1/v3.1.4_TNT_strictConsensus.nwk")
e4_5_conservative = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/5_conservative/round7/v4.7.4_TNT_strictConsensus.nwk")
e4_6_equal = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/6_equal/round4/v5.4.4_TNT_strictConsensus.nwk")
# Create a list of consensus trees 4
e4 = c(e4_1_noPrep, e4_2_noPart, e4_3_max, e4_4_balanced, e4_5_conservative, e4_6_equal)
# Load MPTs
mpt4_1_noPrep = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/1_noPreprocessing/v1.3_PHYG_Domburg_bestUnique.nwk")
mpt4_2_noPart = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/2_noPartitioning/v2.3_PHYG_Domburg_bestUnique.nwk")
mpt4_3_max = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/3_max/v6.3_PHYG_Domburg_bestUnique.nwk")
mpt4_4_balanced = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/4_balanced/round1/v3.1.3_PHYG_Domburg_bestUnique.nwk")
mpt4_5_conservative = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/5_conservative/round7/v4.7.3_PHYG_Domburg_bestUnique.nwk")
mpt4_6_equal = ape::read.tree(file = "../../../Other/IMPORTANT_empirical_data/empirical4_Domburg/6_equal/round4/v5.4.3_PHYG_Domburg_bestUnique.nwk")
# Create a list of MPTs 4
mpt4 = c(mpt4_1_noPrep, mpt4_2_noPart, mpt4_3_max, mpt4_4_balanced, mpt4_5_conservative, mpt4_6_equal)

# Create a list of all trees
all = c(e1, e2, e3, e4)
class(all)

##################################
# PAIRWISE TOPOLOGICAL DISTANCES #
##################################
my_breaks <- seq(0, 1, length.out = 100)

# RF distances in trees 1 (de Sá)
rf1 = TreeDist::RobinsonFoulds(e1, normalize=T)
mat1 <- as.matrix(rf1)
heat1 = pheatmap(mat1, 
         display_numbers = T,     # Shows the values
         number_color = "black", 
         fontsize_number = 10, breaks = my_breaks,
         cluster_rows = F, cluster_cols = F)

# RF distances in trees 2 (Chiari)
rf2 = TreeDist::RobinsonFoulds(e2, normalize=T)
mat2 <- as.matrix(rf2)
heat2 = pheatmap(mat2, 
                 display_numbers = TRUE,     # Shows the values
                 number_color = "black", 
                 fontsize_number = 10,breaks = my_breaks,
                 cluster_rows = F, cluster_cols = F)

# RF distances in trees 3 (Nakamura)
rf3 = TreeDist::RobinsonFoulds(e3, normalize=T)
mat3 <- as.matrix(rf3)
heat3 = pheatmap(mat3, 
                 display_numbers = TRUE,     # Shows the values
                 number_color = "black", 
                 fontsize_number = 10,breaks = my_breaks,
                 cluster_rows = F, cluster_cols = F)

# RF distances in trees 4 (Domburg)
rf4 = TreeDist::RobinsonFoulds(e4, normalize=T)
mat4 <- as.matrix(rf4)
heat4 = pheatmap(mat4, 
                 display_numbers = TRUE,     # Shows the values
                 number_color = "black", 
                 fontsize_number = 10, breaks = my_breaks,
                 cluster_rows = F, cluster_cols = F)

# Heat map
# Convert each one
p1 <- as.ggplot(heat1)
p2 <- as.ggplot(heat2)
p3 <- as.ggplot(heat3)
p4 <- as.ggplot(heat4)
# Now the + and / operators will work!
combined_plot <- (p1 + p2 + p3 + p4) + plot_layout(ncol = 4) +
  plot_annotation(tag_levels = 'a')
# Display it
print(combined_plot)
# Save the plot
ggsave("../empirical_R/heatmaps_RF.pdf", 
       plot = combined_plot, 
       width = 16,      # Width in inches
       height = 4,     # Height in inches
       device = "pdf")


# Load csv with tree costs and rf
data <- read.csv("../empirical_R/empirical_costsXrf.csv")
# Create the new column 'x'
data$x <- log(data$cost_difference)
data$y = log(data$rf_difference)
# General scatter plot
library(ggpubr)
ggscatter(data, x = "x", y = "rf_difference",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue", fill = "lightgray")) +
  stat_cor(method = "pearson", label.x = 3, label.y = 1)   # Automatically adds R and p-value
# Scatter plot with regression line per dataset
ggplot(data, aes(x = cost_difference, y = rf_difference)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkred") +
  theme_minimal() +
  # This creates one plot per 'dataset' category
  facet_wrap(~dataset, scales = "free_x") + 
  labs(title = "RF vs Cost Difference by Dataset",
       x = "Cost Difference",
       y = "RF Difference")

######################################
# PAIRWISE BRANCH LENGTH COMPARISONS #
######################################

###
# Step 1. Use mapBranchLength to map PhyG unique MPTs to PhyG consensus
mpt1_1_noPrep_length = mapBranchLength(e1_1_noPrep, mpt1_1_noPrep)
mpt1_2_noPart_length = mapBranchLength(e1_2_noPart, mpt1_2_noPart)
mpt1_3_max_length = mapBranchLength(e1_3_max, mpt1_3_max)
mpt1_4_balanced_length = mapBranchLength(e1_4_balanced, mpt1_4_balanced)
mpt1_5_conservative_length = mapBranchLength(e1_5_conservative, mpt1_5_conservative)
mpt1_6_equal_length = mapBranchLength(e1_6_equal, mpt1_6_equal)
mpt1_lengths = c(mpt1_1_noPrep_length, mpt1_2_noPart_length, mpt1_3_max_length, mpt1_4_balanced_length, mpt1_5_conservative_length, mpt1_6_equal_length)

mpt2_1_noPrep_length = mapBranchLength(e2_1_noPrep, mpt2_1_noPrep)
mpt2_2_noPart_length = mapBranchLength(e2_2_noPart, mpt2_2_noPart)
mpt2_3_max_length = mapBranchLength(e2_3_max, mpt2_3_max)
mpt2_4_balanced_length = mapBranchLength(e2_4_balanced, mpt2_4_balanced)
mpt2_5_conservative_length = mapBranchLength(e2_5_conservative, mpt2_5_conservative)
mpt2_6_equal_length = mapBranchLength(e2_6_equal, mpt2_6_equal)
mpt2_lengths = c(mpt2_1_noPrep_length, mpt2_2_noPart_length, mpt2_3_max_length, mpt2_4_balanced_length, mpt2_5_conservative_length, mpt2_6_equal_length)

mpt3_1_noPrep_length = mapBranchLength(e3_1_noPrep, mpt3_1_noPrep)
mpt3_2_noPart_length = mapBranchLength(e3_2_noPart, mpt3_2_noPart)
mpt3_3_max_length = mapBranchLength(e3_3_max, mpt3_3_max)
mpt3_4_balanced_length = mapBranchLength(e3_4_balanced, mpt3_4_balanced)
mpt3_5_conservative_length = mapBranchLength(e3_5_conservative, mpt3_5_conservative)
mpt3_6_equal_length = mapBranchLength(e3_6_equal, mpt3_6_equal)
mpt3_lengths = c(mpt3_1_noPrep_length, mpt3_2_noPart_length, mpt3_3_max_length, mpt3_4_balanced_length, mpt3_5_conservative_length, mpt3_6_equal_length)

mpt4_1_noPrep_length = mapBranchLength(e4_1_noPrep, mpt4_1_noPrep)
mpt4_2_noPart_length = mapBranchLength(e4_2_noPart, mpt4_2_noPart)
mpt4_3_max_length = mapBranchLength(e4_3_max, mpt4_3_max)
mpt4_4_balanced_length = mapBranchLength(e4_4_balanced, mpt4_4_balanced)
mpt4_5_conservative_length = mapBranchLength(e4_5_conservative, mpt4_5_conservative)
mpt4_6_equal_length = mapBranchLength(e4_6_equal, mpt4_6_equal)
mpt4_lengths = c(mpt4_1_noPrep_length, mpt4_2_noPart_length, mpt4_3_max_length, mpt4_4_balanced_length, mpt4_5_conservative_length, mpt4_6_equal_length)

###
# Step 2. Use compareBranchLength in each pair of trees

# Dataset 1
# Get the names or indices of the trees
tree_names <- names(mpt1_lengths)
if(is.null(tree_names)) tree_names <- paste0("tree", 1:6)
# Generate all unique pairs (15 pairs for 6 trees)
pairs <- combn(length(mpt1_lengths), 2)
# Use a list to store the results
comparison_results_1 <- list()
for (i in 1:ncol(pairs)) {
  idx1 <- pairs[1, i]
  idx2 <- pairs[2, i]
  # Create a unique name for this pair
  pair_name <- paste0(tree_names[idx1], "_vs_", tree_names[idx2])
  # Run the command
  comparison_results_1[[pair_name]] <- RNODE::compareBranchLength(
    mpt1_lengths[[idx1]], 
    mpt1_lengths[[idx2]], 
    unique = TRUE)}
# Dataset 2
# Get the names or indices of the trees
tree_names <- names(mpt2_lengths)
if(is.null(tree_names)) tree_names <- paste0("tree", 1:6)
# Generate all unique pairs (15 pairs for 6 trees)
pairs <- combn(length(mpt2_lengths), 2)
# Use a list to store the results
comparison_results_2 <- list()
for (i in 1:ncol(pairs)) {
  idx1 <- pairs[1, i]
  idx2 <- pairs[2, i]
  # Create a unique name for this pair
  pair_name <- paste0(tree_names[idx1], "_vs_", tree_names[idx2])
  # Run the command
  comparison_results_2[[pair_name]] <- RNODE::compareBranchLength(
    mpt2_lengths[[idx1]], 
    mpt2_lengths[[idx2]], 
    unique = TRUE)}
# Dataset 3
# Get the names or indices of the trees
tree_names <- names(mpt3_lengths)
if(is.null(tree_names)) tree_names <- paste0("tree", 1:6)
# Generate all unique pairs (15 pairs for 6 trees)
pairs <- combn(length(mpt3_lengths), 2)
# Use a list to store the results
comparison_results_3 <- list()
for (i in 1:ncol(pairs)) {
  idx1 <- pairs[1, i]
  idx2 <- pairs[2, i]
  # Create a unique name for this pair
  pair_name <- paste0(tree_names[idx1], "_vs_", tree_names[idx2])
  # Run the command
  comparison_results_3[[pair_name]] <- RNODE::compareBranchLength(
    mpt3_lengths[[idx1]], 
    mpt3_lengths[[idx2]], 
    unique = TRUE)}
# Dataset 4
# Get the names or indices of the trees
tree_names <- names(mpt4_lengths)
if(is.null(tree_names)) tree_names <- paste0("tree", 1:6)
# Generate all unique pairs (15 pairs for 6 trees)
pairs <- combn(length(mpt4_lengths), 2)
# Use a list to store the results
comparison_results_4 <- list()
for (i in 1:ncol(pairs)) {
  idx1 <- pairs[1, i]
  idx2 <- pairs[2, i]
  # Create a unique name for this pair
  pair_name <- paste0(tree_names[idx1], "_vs_", tree_names[idx2])
  # Run the command
  comparison_results_4[[pair_name]] <- RNODE::compareBranchLength(
    mpt4_lengths[[idx1]], 
    mpt4_lengths[[idx2]], 
    unique = TRUE)}

###
# Step 3. Proportion of branches with different lengths

# Dataset 1
results_vector1 <- sapply(comparison_results_1, function(df) {
  # Remove rows with NA
  df_clean <- na.omit(df)
  # Check if the dataframe is empty after removing NAs to avoid errors
  if (nrow(df_clean) == 0) return(NA)
  # Identify the columns (using dynamic names based on your structure)
  # We assume col 2 is length for tree A and col 4 is length for tree B
  len1 <- df_clean[[2]] 
  len2 <- df_clean[[4]]
  # Perform the calculation
  match_count <- sum((len1 - len2) == 0)
  proportion_different <- 1 - (match_count / nrow(df_clean))
  return(proportion_different)
})
# Dataset 2
results_vector2 <- sapply(comparison_results_2, function(df) {
  # Remove rows with NA
  df_clean <- na.omit(df)
  # Check if the dataframe is empty after removing NAs to avoid errors
  if (nrow(df_clean) == 0) return(NA)
  # Identify the columns (using dynamic names based on your structure)
  # We assume col 2 is length for tree A and col 4 is length for tree B
  len1 <- df_clean[[2]] 
  len2 <- df_clean[[4]]
  # Perform the calculation
  match_count <- sum((len1 - len2) == 0)
  proportion_different <- 1 - (match_count / nrow(df_clean))
  return(proportion_different)
})
# Dataset 3
results_vector3 <- sapply(comparison_results_3, function(df) {
  # Remove rows with NA
  df_clean <- na.omit(df)
  # Check if the dataframe is empty after removing NAs to avoid errors
  if (nrow(df_clean) == 0) return(NA)
  # Identify the columns (using dynamic names based on your structure)
  # We assume col 2 is length for tree A and col 4 is length for tree B
  len1 <- df_clean[[2]] 
  len2 <- df_clean[[4]]
  # Perform the calculation
  match_count <- sum((len1 - len2) == 0)
  proportion_different <- 1 - (match_count / nrow(df_clean))
  return(proportion_different)
})
# Dataset 4
results_vector4 <- sapply(comparison_results_4, function(df) {
  # Remove rows with NA
  df_clean <- na.omit(df)
  # Check if the dataframe is empty after removing NAs to avoid errors
  if (nrow(df_clean) == 0) return(NA)
  # Identify the columns (using dynamic names based on your structure)
  # We assume col 2 is length for tree A and col 4 is length for tree B
  len1 <- df_clean[[2]] 
  len2 <- df_clean[[4]]
  # Perform the calculation
  match_count <- sum((len1 - len2) == 0)
  proportion_different <- 1 - (match_count / nrow(df_clean))
  return(proportion_different)
})

###
# Step 4. Create a matrix with proportions

# Define a function to transform the vector into a matrix
vector_to_matrix <- function(vec, tree_names = paste0("tree", 1:6)) {
  n <- length(tree_names)
  # Initialize empty matrix
  m <- matrix(0, nrow = n, ncol = n, dimnames = list(tree_names, tree_names))
  # Fill the matrix
  for (name in names(vec)) {
    parts <- strsplit(name, "_vs_")[[1]]
    t1 <- parts[1]
    t2 <- parts[2]
    val <- vec[name]
    m[t1, t2] <- val
    m[t2, t1] <- val
  }
  return(m)
}
# Run the function for each of your objects
mat1_branch <- vector_to_matrix(results_vector1)
mat2_branch <- vector_to_matrix(results_vector2)
mat3_branch <- vector_to_matrix(results_vector3)
mat4_branch <- vector_to_matrix(results_vector4)
# Delete the string tree
rownames(mat1_branch) <- gsub("tree", "", rownames(mat1_branch))
colnames(mat1_branch) <- gsub("tree", "", colnames(mat1_branch))
rownames(mat2_branch) <- gsub("tree", "", rownames(mat2_branch))
colnames(mat2_branch) <- gsub("tree", "", colnames(mat2_branch))
rownames(mat3_branch) <- gsub("tree", "", rownames(mat3_branch))
colnames(mat3_branch) <- gsub("tree", "", colnames(mat3_branch))
rownames(mat4_branch) <- gsub("tree", "", rownames(mat4_branch))
colnames(mat4_branch) <- gsub("tree", "", colnames(mat4_branch))

###
# Step 5. Plot the heatmap

#Dataset 1
labels1 = matrix(sprintf("%.1f%%", mat1_branch * 100),
                nrow = nrow(mat1_branch))
heat1_branch = pheatmap(mat1_branch, 
                 display_numbers = labels1,     # Shows the values
                 number_color = "black", 
                 fontsize_number = 10,breaks = my_breaks,
                 cluster_rows = F, cluster_cols = F)
#Dataset 2
labels2 = matrix(sprintf("%.1f%%", mat2_branch * 100),
                nrow = nrow(mat2_branch))
heat2_branch = pheatmap(mat2_branch, 
                        display_numbers = labels2,     # Shows the values
                        number_color = "black", 
                        fontsize_number = 10,breaks = my_breaks,
                        cluster_rows = F, cluster_cols = F)
#Dataset 3
labels3 = matrix(sprintf("%.1f%%", mat3_branch * 100),
                 nrow = nrow(mat3_branch))
heat3_branch = pheatmap(mat3_branch, 
                        display_numbers = labels3,     # Shows the values
                        number_color = "black", 
                        fontsize_number = 10,breaks = my_breaks,
                        cluster_rows = F, cluster_cols = F)
#Dataset 4
labels4 = matrix(sprintf("%.1f%%", mat4_branch * 100),
                 nrow = nrow(mat4_branch))
heat4_branch = pheatmap(mat4_branch, 
                        display_numbers = labels4,     # Shows the values
                        number_color = "black", 
                        fontsize_number = 10,breaks = my_breaks,
                        cluster_rows = F, cluster_cols = F)
# Heat map
# Convert each one
p1_branch <- as.ggplot(heat1_branch)
p2_branch <- as.ggplot(heat2_branch)
p3_branch <- as.ggplot(heat3_branch)
p4_branch <- as.ggplot(heat4_branch)
# Now the + and / operators will work!
combined_plot_branch <- (p1_branch + p2_branch + p3_branch + p4_branch) + 
  plot_layout(ncol = 4) +
  plot_annotation(tag_levels = 'a')
# Display it
print(combined_plot_branch)
# Save the plot
ggsave("../empirical_R/heatmaps_branchLengths.pdf", 
       plot = combined_plot_branch, 
       width = 16,      # Width in inches
       height = 4,     # Height in inches
       device = "pdf")

##############################
# PLOT TWO TREES AS EXAMPLES #
##############################

# Plot two trees with unique nodes as example
RNODE::uniqueNodes(mpt2_1_noPrep_length, mpt2_6_equal_length, 
                   plotTrees = T, node.numbers = F, sup=F, tree.fsize=1, tree.height=8,
                   output.tree="../empirical_R/fig_Chiari_1x6.pdf")

# Plot two trees without unique nodes as example (add values of branch lengths in Photoshop)
RNODE::uniqueNodes(mpt2_3_max_length, mpt2_6_equal_length, 
                   plotTrees = T, node.numbers = F, sup=F, tree.fsize=1, tree.height=8,
                   output.tree="../empirical_R/fig_Chiari_3x6.pdf")
a=RNODE::compareBranchLength(mpt2_3_max_length, mpt2_6_equal_length, composition=T)

#!/usr/bin/env Rscript

#############
# PRE-AMBLE #
#############

# Set the working directory 
setwd("/Users/labanfibios/Desktop/Doutorado/Project/B3_PrepDyn/GitHub/test_data/simulations_v1_R/")

# Load required packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)
if (!require(GGally)) install.packages("GGally", dependencies = TRUE)
if (!require(viridis)) install.packages("viridis")
if (!require(readxl)) install.packages("readxl", dependencies = TRUE)

#############
# LOAD DATA #
#############

# Read the data
df <- read_excel("simulations_complexity.xlsx", sheet = "new")

# Inspect the data structure
str(df)

#################
# PREPROCESSING #
#################

# Run model
model <- lm(log(prepDyn_CPU_time) ~ log(n_leaves) + log(n_columns) + log(n_partitions), data = df)
summary(model)

# Plot: Time vs Partitions, grouped by (n_leaves, n_columns)
df$group <- interaction(df$n_leaves, df$n_columns, drop = TRUE)
p <- ggplot(df, aes(x = n_partitions, y = prepDyn_CPU_time, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "No. partitions",
       y = "CPU time (s) for preprocessing",
       color = "No. leaves × No. characters") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
ggsave("simulations_complexityPartitions_prepDyn.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

# Plot: Time vs Leaves, grouped by (n_leaves, n_columns)
df$group <- interaction(df$n_columns, df$n_partitions, drop = TRUE)
p <- ggplot(df, aes(x = n_leaves, y = prepDyn_CPU_time, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "No. leaves",
       y = "CPU time (s) for preprocessing",
       color = "No. characters × No. partitions") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
ggsave("simulations_complexityLeaves_prepDyn.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

# Plot: Time vs Characters, grouped by (n_leaves, n_partitions)
df$group <- interaction(df$n_leaves, df$n_partitions, drop = TRUE)
p <- ggplot(df, aes(x = n_columns, y = prepDyn_CPU_time, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "No. characters",
       y = "CPU time (s) for preprocessing",
       color = "No. leaves × No. partitions") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
ggsave("simulations_complexityCharacters_prepDyn.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

sink()

#########################
# PHYLOGENETIC ANALYSES #
#########################

# Run model
model <- lm(log(cost) ~ log(n_leaves) + log(n_columns) + log(n_partitions), data = df)
summary(model)

# Run model
model <- lm(log(swap_CPU_time) ~ log(n_leaves) + log(n_columns) + log(n_partitions), data = df)
summary(model)

# Plot: Cost vs Partitions, grouped by (n_leaves, n_columns)
df$group <- interaction(df$n_leaves, df$n_columns, drop = TRUE)
p <- ggplot(df, aes(x = n_partitions, y = cost, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "No. partitions",
       y = "Parsimony score",
       color = "No. leaves × No. characters") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
ggsave("simulations_empiricalCostPartitions_prepDyn.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

# Plot: Cost vs Partitions, grouped by (n_leaves, n_columns)
df$group <- interaction(df$n_leaves, df$n_columns, drop = TRUE)
p <- ggplot(df, aes(x = n_partitions, y = swap_CPU_time, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "No. partitions",
       y = "CPU time (s) for phylogenetic analyses",
       color = "No. leaves × No. characters") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
ggsave("simulations_empiricalTimePartitions_prepDyn.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

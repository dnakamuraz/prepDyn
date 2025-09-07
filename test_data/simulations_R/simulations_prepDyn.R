#!/usr/bin/env Rscript

# Load required packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)
if (!require(GGally)) install.packages("GGally", dependencies = TRUE)
if (!require(viridis)) install.packages("viridis")

# Load libraries
library(ggplot2)
library(GGally)
library(viridis)

# Read the data
df <- read.csv("simulations_complexity.csv")

# Inspect the data structure
str(df)

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
       y = "CPU time (s)",
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
       y = "CPU time (s)",
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
       y = "CPU time (s)",
       color = "No. leaves × No. partitions") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
ggsave("simulations_complexityCharacters_prepDyn.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

sink()

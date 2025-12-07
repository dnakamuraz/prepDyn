#!/usr/bin/env Rscript

#############
# PRE-AMBLE #
#############

# Set the working directory 
setwd("/Users/labanfibios/Desktop/Doutorado/Project/B3_PrepDyn/GitHub/test_data/simulations_v1_R/")

# Load required packages
if (!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)
if (!require(GGally)) install.packages("GGally", dependencies = TRUE)
if (!require(mgcv)) install.packages("mgcv", dependencies = TRUE)
if (!require(viridis)) install.packages("viridis")
if (!require(readxl)) install.packages("readxl", dependencies = TRUE)

#############
# LOAD DATA #
#############

# Read the data
df <- read_excel("simulations_complexity.xlsx", sheet = "with80t")

# Inspect the data structure
str(df)

#################
# PREPROCESSING #
#################

# Run model
model <- lm(log2(prepDyn_CPU_time) ~ log2(n_leaves) + log2(n_columns) + log2(n_partitions), data = df)
summary(model)

# Plot: Time vs Partitions, grouped by (n_leaves, n_columns)
df$group <- interaction(df$n_leaves, df$n_columns, drop = TRUE)
p <- ggplot(df, aes(x = n_partitions, y = prepDyn_CPU_time, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "\nNo. partitions",
       y = "CPU time (s) for preprocessing\n",
       color = "No. leaves × No. nucleotides") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
p
ggsave("simulations_noPartitionsXtimePreprocessing.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

# Plot: Time vs Leaves, grouped by (n_leaves, n_columns)
df$group <- interaction(df$n_columns, df$n_partitions, drop = TRUE)
p <- ggplot(df, aes(x = n_leaves, y = prepDyn_CPU_time, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "\nNo. leaves",
       y = "CPU time (s) for preprocessing\n",
       color = "No. nucleotides × No. partitions") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
p
ggsave("simulations_noLeavesXtimePreprocessing.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

# Plot: Time vs Characters, grouped by (n_leaves, n_partitions)
df$group <- interaction(df$n_leaves, df$n_partitions, drop = TRUE)
p <- ggplot(df, aes(x = n_columns, y = prepDyn_CPU_time, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "\nNo. nucleotides",
       y = "CPU time (s) for preprocessing\n",
       color = "No. leaves × No. partitions") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
p
ggsave("simulations_noCharactersXtimePreprocessing.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

sink()

###################################
# PHYLOGENETIC ANALYSES: DO COSTS #
###################################

# Cost
model <- lm(log2(cost) ~ log2(n_leaves) + log2(n_columns) + log2(n_partitions), data = df)
summary(model)

# Plot: Cost vs Partitions, grouped by (n_leaves, n_columns)
df$group <- interaction(df$n_leaves, df$n_columns, drop = TRUE)
p <- ggplot(df, aes(x = n_partitions, y = cost, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "\nNo. partitions",
       y = "Parsimony score\n",
       color = "No. leaves × No. nucleotides") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
p
ggsave("simulations_noPartitionsXcostDO.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

##################################
# PHYLOENETIC ANALYSES: DO TIME #
##################################

# LINEAR MODELS
# Runtime ~ leaves * characters * partitions 
model1 <- lm(log2(swap_CPU_time) ~ log2(n_leaves) * log2(n_columns) * log2(n_partitions), data = df)
summary(model1)
# Runtime ~ leaves + characters + partitions + leaves:partitions + characters:partitions
model2 <- lm(log2(swap_CPU_time) ~ log2(n_leaves) + log2(n_columns) + log2(n_partitions) + log2(n_leaves):log2(n_partitions) + log2(n_columns):log2(n_partitions), data = df)
summary(model2)
# Runtime ~ leaves + characters + partitions + leaves:partitions
model3 <- lm(log2(swap_CPU_time) ~ log2(n_leaves) + log2(n_columns) + log2(n_partitions) + log2(n_leaves):log2(n_partitions), data = df)
summary(model3)
# Runtime ~ leaves + characters + partitions + leaves:characters
model4 <- lm(log2(swap_CPU_time) ~ log2(n_leaves) + log2(n_columns) + log2(n_partitions) + log2(n_columns):log2(n_partitions), data = df)
summary(model4)
# Runtime ~ leaves + characters + partitions
model5 <- lm(log2(swap_CPU_time) ~ log2(n_leaves) + log2(n_columns) + log2(n_partitions), data = df)
summary(model5)

# NONLINEAR MODELS
# Runtime ~ leaves * characters * partitions 
model6 = gam(log2(swap_CPU_time) ~ s(log2(n_leaves),k=4) + 
               s(log2(n_columns),k=3) + 
               s(log2(n_partitions),k=8) +
               te(log2(n_leaves),log2(n_partitions), k=c(4,8)) +
               te(log2(n_columns),log2(n_partitions), k=c(3,8)) +
               te(log2(n_columns),log2(n_partitions),log2(n_leaves), k=c(3,8,4)), 
             data=df)
summary(model6)
# Runtime ~ leaves + characters + partitions + leaves:partitions + characters:partitions
model7 = gam(log2(swap_CPU_time) ~ s(log2(n_leaves),k=4) + 
               s(log2(n_columns),k=3) + 
               s(log2(n_partitions),k=8) +
               te(log2(n_leaves),log2(n_partitions), k=c(4,8)) +
               te(log2(n_columns),log2(n_partitions), k=c(3,8)), 
             data=df)
summary(model7)
# Runtime ~ leaves + characters + partitions + leaves:partitions
model8 = gam(log2(swap_CPU_time) ~ s(log2(n_leaves),k=4) + 
               s(log2(n_columns),k=3) + 
               s(log2(n_partitions),k=8) +
               te(log2(n_leaves),log2(n_partitions), k=c(4,8)), 
             data=df)
summary(model8)
# Runtime ~ leaves + characters + partitions + characters:partitions
model9 = gam(log2(swap_CPU_time) ~ s(log2(n_leaves),k=4) + 
               s(log2(n_columns),k=3) + 
               s(log2(n_partitions),k=8) +
               te(log2(n_columns),log2(n_partitions), k=c(3,8)), 
             data=df)
summary(model9)
# Runtime ~ leaves + characters + partitions
model10 = gam(log2(swap_CPU_time) ~ s(log2(n_leaves),k=4) + 
               s(log2(n_columns),k=3) + 
               s(log2(n_partitions),k=8), 
             data=df)
summary(model10)

# MODEL SELECTION
# AIC of linear models
AIC(model1, model2, model3, model4, model5)
# AIC of linear models
AIC(model6, model7, model8, model9, model10)
# AIC of all models
AIC(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)

# VISUALIZATION
# Plot: time vs Partitions, grouped by (n_leaves, n_columns)
df$group <- interaction(df$n_leaves, df$n_columns, drop = TRUE)
p <- ggplot(df, aes(x = n_partitions, y = swap_CPU_time, color = group)) +
  geom_line() +
  geom_point() +
  labs(title = "",
       x = "\nNo. partitions",
       y = "CPU time (s) for direct optimization\n",
       color = "No. leaves × No. nucleotides") +
  theme_classic() +
  scale_color_viridis_d(option = "D")  # discrete viridis palette
p
ggsave("simulations_noPartitionsXtimeDO.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

# 2D plot
grid <- expand.grid(
  n_leaves = median(df$n_leaves),
  n_columns = seq(min(df$n_columns), max(df$n_columns), length=50),
  n_partitions = seq(min(df$n_partitions), max(df$n_partitions), length=50)
)
# Predict CPU time
grid$pred <- predict(model4, newdata=grid)
grid$pred_time <- 2^grid$pred  # back to CPU time
# Contour plot
p=ggplot(grid, aes(x = n_columns, y = n_partitions, z = pred_time)) +
  geom_contour_filled() +
  scale_x_log10() + scale_y_log10() +
  labs(x="\nNo. nucleotides", y="No. partitions\n", fill="CPU time") +
  theme_minimal()
p
ggsave("simulations_interactionPartitionsNucleotidesXtimeDO.jpg", p, width = 8, height = 5, dpi = 300, units = "in")

# 3D plot
grid <- expand.grid(
  n_leaves = median(df$n_leaves),                # hold leaves constant
  n_columns = seq(min(df$n_columns), max(df$n_columns), length.out = 100),
  n_partitions = seq(min(df$n_partitions), max(df$n_partitions), length.out = 100)
)
grid$pred <- predict(model4, newdata = grid)
grid$CPU_time <- 2^grid$pred
z_matrix <- matrix(grid$CPU_time, nrow=100, ncol=100)
library(plotly)
plot_ly(
  x = seq(min(df$n_columns), max(df$n_columns), length.out=100),
  y = seq(min(df$n_partitions), max(df$n_partitions), length.out=100),
  z = z_matrix,
  type = "surface"
) %>%
  layout(
    scene = list(
      xaxis = list(title="Number of columns", type="log"),
      yaxis = list(title="Number of partitions", type="log"),
      zaxis = list(title="CPU time")
    )
  )

####################################################################################
# PHYLOGENETIC ANALYSES: DO TIME USING A SUBSAMPLE (40 LEAVES + 10000 NUCLEOTIDES) #
####################################################################################

df_subset <- subset(df, n_leaves == 40 & n_columns == 10000)

# LINEAR MODEL
model_lm <- lm(log2(swap_CPU_time) ~ log2(n_partitions), data = df_subset)
summary(model_lm)
# NONLINEAR MODEL
#Gam
model_gam = gam(log2(swap_CPU_time) ~ s(log2(n_partitions),k=8), 
              data=df_subset)
summary(model_gam)
#Polynomial
model_poly1 <- lm(log2(swap_CPU_time) ~ poly(log2(n_partitions), 1), data = df_subset)
model_poly2 <- lm(log2(swap_CPU_time) ~ poly(log2(n_partitions), 2), data = df_subset)
model_poly3 <- lm(log2(swap_CPU_time) ~ poly(log2(n_partitions), 3), data = df_subset)
model_poly4 <- lm(log2(swap_CPU_time) ~ poly(log2(n_partitions), 4), data = df_subset)
AIC(model_poly1, model_poly2, model_poly3, model_poly4)
summary(model_poly4)

# AIC
AIC(model_lm, model_gam, model_poly4)

# PLOT GAM
library(mgcv)
# Base plot
plot(model_gam, shade=TRUE, main="Effect of n_partitions on CPU time")
# Or ggplot with predicted values
newdata <- data.frame(n_partitions = seq(min(df_subset$n_partitions),
                                         max(df_subset$n_partitions),
                                         length.out=100))
newdata$pred <- predict(model_gam, newdata=newdata, type="response")
newdata$CPU_time <- 2^newdata$pred  # back to original scale
ggplot(newdata, aes(x=n_partitions, y=CPU_time)) +
  geom_line() +
  labs(x="Number of partitions", y="CPU time") +
  theme_minimal()

# PLOT POLYNOMIAL
# Create a sequence of partition values
newdata <- data.frame(n_partitions = seq(min(df_subset$n_partitions),
                                         max(df_subset$n_partitions),
                                         length.out=100))
# Predict log2(CPU_time)
newdata$pred <- predict(model_poly4, newdata=newdata)
newdata$CPU_time <- 2^newdata$pred  # back to original scale
# Plot
library(ggplot2)
ggplot(newdata, aes(x=n_partitions, y=CPU_time)) +
  geom_line(color="blue", size=1) +
  geom_point(data=df_subset, aes(x=n_partitions, y=swap_CPU_time), alpha=0.5) +
  labs(x="Number of partitions", y="CPU time", 
       title="Polynomial regression of CPU time vs partitions") +
  theme_minimal()

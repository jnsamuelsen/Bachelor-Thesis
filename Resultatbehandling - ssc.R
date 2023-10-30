# Clearing global environment
rm(list = ls())

# -------------------------------------------------------------------------------------------------------------------
# Loading libraries

library(tidyverse)
library(ggplot2)
library(GGally)
library(dplyr)

# -------------------------------------------------------------------------------------------------------------------
# Loading cell subset data into R
data <- read.csv("C:/Users/jeppe/Desktop/Bachelorprojekt/CIBERSORT/Results/ssc_results.csv")

# -------------------------------------------------------------------------------------------------------------------
# Data wrangling

# Converting to tibble
tibble_data <- as.tibble(data)

# Showing tibble
tibble_data

# Creating a vector that maps "h" and "s" mixtures to a healthy and ill status_marker
status_marker <- c()
for (i in tibble_data$Mixture) {
  if (str_detect(i, "^h")) {
    status_marker <- c(status_marker, "healthy")
  } 
  else if (str_detect(i, "^s")) {
    status_marker <- c(status_marker, "ill")
  }
}

# We know want to add a column label that indicates the status for the patient (either healthy or ill)
tibble_data <- tibble_data %>% 
  mutate(status = status_marker)

# -------------------------------------------------------------------------------------------------------------------
### Boxplots ###

# Boxplot of CD8 T-cells, facetted on status
ggplot(data = tibble_data, mapping = aes(x = "", y = T.cells.CD8)) +
  geom_boxplot() +
  facet_wrap(~status)

# Boxplot of memory B-cells, facetted on status
ggplot(data = tibble_data, mapping = aes(x = "", y = B.cells.memory)) +
  geom_boxplot() +
  facet_wrap(~status)

# Boxplot of CD4 naive T-cells, facetted on status
ggplot(data = tibble_data, mapping = aes(x = "", y = T.cells.CD4.naive)) +
  geom_boxplot() +
  facet_wrap(~status)

# Boxplot of memory activated CD4 memory T-cells, facetted on status
ggplot(data = tibble_data, mapping = aes(x = "", y = T.cells.CD4.memory.activated)) +
  geom_boxplot() +
  facet_wrap(~status)

# Boxplot of resting NK cells, facetted on status
ggplot(data = tibble_data, mapping = aes(x = "", y = NK.cells.resting)) +
  geom_boxplot() +
  facet_wrap(~status)

# Boxplot of activated NK cells, facetted on status
ggplot(data = tibble_data, mapping = aes(x = "", y = NK.cells.activated)) +
  geom_boxplot() +
  facet_wrap(~status)

# Boxplot of monocytes, facetted on status
ggplot(data = tibble_data, mapping = aes(x = "", y = Monocytes)) +
  geom_boxplot() +
  facet_wrap(~status)

# Boxplot of activated dendritic cells, facetted on status
ggplot(data = tibble_data, mapping = aes(x = "", y = Dendritic.cells.activated)) +
  geom_boxplot() +
  facet_wrap(~status)

# Boxplot of neutrophils, facetted on status
ggplot(data = tibble_data, mapping = aes(x = "", y = Neutrophils)) +
  geom_boxplot() +
  facet_wrap(~status)

# -------------------------------------------------------------------------------------------------------------------
### Stacked bar plots ###


# Splitting the data into two tibbles one for ill and one for healthy individuals
tibble_ill <- tibble_data %>% filter(status == "ill")
tibble_healthy <-  tibble_data %>% filter(status == "healthy")


# Now creating vectors which holds mean values of fractions for the different immune cells for respectively healthy and ill individuals
healthy_means <- c()
ill_means <- c()

for (col_ill in tibble_ill[2:23]) {
  ill_means <- c(ill_means, mean(col_ill))
}

for (col_healthy in tibble_healthy[2:23]) {
  healthy_means <- c(healthy_means, mean(col_healthy))
}

# We now create a new tibble to play with cell_type in first column, mean_cell fraction in 2nd column and health status in 3rd column
cell_tibble <- tibble(cell_type = c(colnames(tibble_data)[2:23], colnames(tibble_data)[2:23]), mean_fraction = c(ill_means, healthy_means), health_status = c(rep("ill",22),rep("healthy",22)))

# Creating where the neutrophils are filtered away
cell_tibble_no_neutrophils <- cell_tibble %>% 
  filter(mean_fraction < 0.5)

# Plotting healthy vs ill stacked bar plot (each type mean fraction stacked on top of each other)
ggplot(data = cell_tibble, aes(x = health_status, y = mean_fraction, fill = cell_type)) +
  geom_col()

# Plotting healthy vs ill stacked bar plot (each type mean fraction stacked on top of each other) - without neutrophils
ggplot(data = cell_tibble_no_neutrophils, aes(x = health_status, y = mean_fraction, fill = cell_type)) +
  geom_col()

# Plotting each cell type mean fraction bar plot (with healthy vs ill fractions stacked on top of each other) - without neutrophils
ggplot(data = cell_tibble_no_neutrophils, aes(x = cell_type, y = mean_fraction, fill = health_status)) +
  geom_col() +
  coord_flip()

# -------------------------------------------------------------------------------------------------------------------
### Mann-Whitney U test (Wilcoxon test) ###

# We will start by doing the test for the T-cells CD8
tcd8_wilcoxon_test <- wilcox.test(tibble_ill$T.cells.CD8, tibble_healthy$T.cells.CD8)
tcd8_wilcoxon_test
# At a significance level of 5 %, we can conclude that cells in healthy individuals are significantly differential than in ssc patients
group_by(tibble_data, status) %>% 
  summarise(
    count = n(),
    median = median(T.cells.CD8, na.rm = TRUE),
    IQR = IQR(T.cells.CD8, na.rm = TRUE)
  )
# Based on the summary statistics we realise that T-cells CD8 are significantly higher expressed in healthy individuals than in ssc patients

# Now doing the test for B-cells memory
bmemory_wilcoxon_test <- wilcox.test(tibble_ill$B.cells.memory, tibble_healthy$B.cells.memory)
bmemory_wilcoxon_test
# At a significance level of 5 %, we can conclude that cells in healthy individuals aren't significantly different than in ssc patients

# Now doing the test for T-cells CD4 naive
tcd4naive_wilcoxon_test <- wilcox.test(tibble_ill$T.cells.CD4.naive, tibble_healthy$T.cells.CD4.naive)
tcd4naive_wilcoxon_test
# At a significance level of 5 %, we can conclude that cells in healthy individuals are significantly different than in ssc patients
group_by(tibble_data, status) %>% 
  summarise(
    count = n(),
    median = median(T.cells.CD4.naive, na.rm = TRUE),
    IQR = IQR(T.cells.CD4.naive, na.rm = TRUE)
  )
# Based on the summary statistics we realise that T-cells CD4 naive are significantly higher expressed in healthy individuals than in ssc patients

# Now doing the test for Macrophages M0
tcellsgammadelta_wilcoxon_test <- wilcox.test(tibble_ill$T.cells.gamma.delta, tibble_healthy$T.cells.gamma.delta)
tcellsgammadelta_wilcoxon_test
# At a significance level of 5 %, we can conclude that cells in healthy individuals aren't significantly different than in ssc patients



# -------------------------------------------------------------------------------------------------------------------
### Heatmap ###


# -------------------------------------------------------------------------------------------------------------------
# Correlations between individual cell types?

# This can be done like the matrix plot in ML
# Correlation between cd8 t-cells and cd4 naive t-cells, colored by healthy/ill?
ggplot(data = tibble_data, mapping = aes(x = T.cells.CD4.naive, y = T.cells.CD8, color = status)) +
  geom_point()

# Histogram
ggplot(data = tibble_data, mapping = aes(x = T.cells.CD8)) +
  geom_histogram() +
  facet_wrap(~status)

# -------------------------------------------------------------------------------------------------------------------
# Statistical tests of independence

# Create boxplots for every cell type, "healthy" vs "ill", + outlier detection

# Bar plot of healthy vs ill for every cell type (mean value of e.g. healthy status and T cells CD8)

# Correlation between individual cell types

# Stacked bar plot

# Heat map?

# Other plots?

# Mann-Whitney U test of statistical significance
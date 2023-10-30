# Clearing global environment
rm(list = ls())

# Loading libraries
library(tidyverse)
library(broom)
library(ggplot2)

# Loading cell subset data into R
data <- read.csv("C:/Users/jeppe/Desktop/Bachelorprojekt/CIBERSORT/Results/ssc_results.csv")

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

# Selecting only the numerical columns that we need for the pca
pca_tibble <- tibble_data %>% 
  select(-Mixture, -P.value, -Correlation, -RMSE, -status)

# Calculating principal components
ssc_pca <- pca_tibble %>% 
  prcomp(center = TRUE, scale. = FALSE)

# Using broom to tidy the pca data
ssc_pca %>% tidy("pcs")

# Plotting the variance explained by the principal components
ssc_pca %>%
  tidy("pcs") %>% 
  ggplot(mapping = aes(x = PC, y = percent)) +
  geom_col() +
  theme_bw()

ssc_pca %>% tidy("samples")

# Using broom to tidy and augment ssc data
ssc_pca_aug <- ssc_pca %>% augment(tibble_data)

# Plotting the principal components
ssc_pca_aug %>% 
  ggplot(mapping = aes(x = .fittedPC1, y = .fittedPC2, colour = status)) +
  geom_point()

# Clearing global environment
rm(list = ls())

# -------------------------------------------------------------------------------------------------------------------
# Loading libraries

library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(xtable)
library(tidyr)

# Loading data into R
load("C:/Users/jeppe/Desktop/Bachelorprojekt/Data/MS data/data_updated.RData")

# Transposing for tidy data and converting to tibbles
xCell_matrix_adj <- t(xCell_matrix_adj)

xCell_tibble_adj <- xCell_matrix_adj %>% 
  as_tibble()

xCell_matrix <- t(xCell_matrix)

xCell_tibble <- xCell_matrix %>% 
  as_tibble()

xCell_matrix_age_adj <- t(xCell_matrix_age_adj)

xCell_tibble_age_adj <- xCell_matrix_age_adj %>% 
  as_tibble()

# Renaming pheno columns
pheno <- pheno %>% 
  rename(Status = R.RRMS..S.HC..G.Gilenya.behandlet.,
         Sex = sex.male1,
         Past_treatment = Tidligere.Beh,
         Stability = X0Stabil_1ustabil)

# Adding the filename column to the tibbles
xCell_tibble_adj <- xCell_tibble_adj %>%
  mutate(Filename = rownames(xCell_matrix_adj))

xCell_tibble <- xCell_tibble %>% 
  mutate(Filename = rownames(xCell_matrix))

xCell_tibble_age_adj <- xCell_tibble_age_adj %>%
  mutate(Filename = rownames(xCell_matrix_age_adj))

# Joining the pheno data to the xCell tibbles
xCell_tibble_adj <- full_join(xCell_tibble_adj, pheno, by = "Filename")
xCell_tibble <- full_join(xCell_tibble, pheno, by = "Filename")
xCell_tibble_age_adj <- full_join(xCell_tibble_age_adj, pheno, by = "Filename")

# -------------------------------------------------------------------------------------------------------------------
## Investigating distributions of enrichment scores

## RRMS vs HC density plots
# CD8+ Tcm
Distri_R_S_CD8plus_Tcm <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = `CD8+ Tcm`, fill = Status)) +
  geom_density(alpha = .5) +
  labs(title = "CD8+ Tcm cell enrichment distributions", x = "CD8+ Tcm xCell enrichment score", y = "Density", fill = "Health status") +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 16))

# CD8+ T-cells
Distri_R_S_CD8plus_Tcells <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = `CD8+ T-cells`, fill = Status)) +
  geom_density(alpha = .5) +
  labs(title = "CD8+ T-cell enrichment distributions", x = "CD8+ T-cells xCell enrichment score", y = "Density", fill = "Health status") +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 16))

# Th1 cells
Distri_R_S_Th1_cells <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = `Th1 cells`, fill = Status)) +
  geom_density(alpha = .5) +
  labs(title = "Th1 cell enrichment distributions", x = "Th1 cells xCell enrichment score", y = "Density", fill = "Health status") +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 16))

# Th2 cells
Distri_R_S_Th2_cells <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = `Th2 cells`, fill = Status)) +
  geom_density(alpha = .5) +
  labs(title = "Th2 cell enrichment distributions", x = "Th2 cells xCell enrichment score", y = "Density", fill = "Health status") +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 16))

## Active vs inactive on fingolimod

# Renaming Stability column values to 0 ~ stabile, 1 ~ unstabile
xCell_tibble_adj <- xCell_tibble_adj %>% 
  mutate(Stability = case_when(Stability == 0 ~ "Stabile",
                               Stability == 1 ~ "Unstabile"))

# CD8+ Tcm
Distri_G_CD8plus_Tcm <- xCell_tibble_adj %>% 
  filter(Status == "G") %>% 
  ggplot(mapping = aes(x = `CD8+ Tcm`, fill = Stability)) +
  geom_density(alpha = .5) +
  labs(title = "CD8+ Tcm cell enrichment distributions", x = "CD8+ Tcm xCell enrichment score", y = "Density", fill = "MS patient stability") +
  scale_fill_brewer(palette = "Dark2") +
  theme(text = element_text(size = 16))

# CD8+ T-cells
Distri_G_CD8plus_Tcells <- xCell_tibble_adj %>% 
  filter(Status == "G") %>% 
  ggplot(mapping = aes(x = `CD8+ T-cells`, fill = Stability)) +
  geom_density(alpha = .5) +
  labs(title = "CD8+ T-cell enrichment distributions", x = "CD8+ T-cells xCell enrichment score", y = "Density", fill = "MS patient stability") +
  scale_fill_brewer(palette = "Dark2") +
  theme(text = element_text(size = 16))

# Th1 cells
Distri_G_Th1_cells <- xCell_tibble_adj %>% 
  filter(Status == "G") %>% 
  ggplot(mapping = aes(x = `Th1 cells`, fill = Stability)) +
  geom_density(alpha = .5) +
  labs(title = "Th1 cell enrichment distributions", x = "Th1 cells xCell enrichment score", y = "Density", fill = "MS patient stability") +
  scale_fill_brewer(palette = "Dark2") +
  theme(text = element_text(size = 16))

# Th2 cells
Distri_G_Th2_cells <- xCell_tibble_adj %>% 
  filter(Status == "G") %>% 
  ggplot(mapping = aes(x = `Th2 cells`, fill = Stability)) +
  geom_density(alpha = .5) +
  labs(title = "Th2 cell enrichment distributions", x = "Th2 cells xCell enrichment score", y = "Density", fill = "MS patient stability") +
  scale_fill_brewer(palette = "Dark2") +
  theme(text = element_text(size = 16))

## RRMS females vs RRMS males

# Renaming Sex column values to 0 ~ Female, 1 ~ Male
xCell_tibble_age_adj <- xCell_tibble_age_adj %>% 
  mutate(Sex = case_when(Sex == 0 ~ "Female",
                         Sex == 1 ~ "Male"))

# CD8+ Tcm
Distri_R_S_gender_CD8plus_Tcm <- xCell_tibble_age_adj %>% 
  filter(Status == "R") %>% 
  ggplot(mapping = aes(x = `CD8+ Tcm`, fill = Sex)) +
  geom_density(alpha = .5) +
  labs(title = "CD8+ Tcm cell enrichment distributions", x = "CD8+ Tcm xCell enrichment score", y = "Density", fill = "Gender") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 16))

# CD8+ T-cells
Distri_R_S_gender_CD8plus_Tcells <- xCell_tibble_age_adj %>% 
  filter(Status == "R") %>% 
  ggplot(mapping = aes(x = `CD8+ T-cells`, fill = Sex)) +
  geom_density(alpha = .5) +
  labs(title = "CD8+ T-cell enrichment distributions", x = "CD8+ T-cells xCell enrichment score", y = "Density", fill = "Gender") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 16))

# Th1 cells
Distri_R_S_gender_Th1_cells <- xCell_tibble_age_adj %>% 
  filter(Status == "R") %>% 
  ggplot(mapping = aes(x = `Th1 cells`, fill = Sex)) +
  geom_density(alpha = .5) +
  labs(title = "Th1 cell enrichment distributions", x = "Th1 cells xCell enrichment score", y = "Density", fill = "Gender") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 16))

# Th2 cells
Distri_R_S_gender_Th2_cells <- xCell_tibble_age_adj %>% 
  filter(Status == "R") %>% 
  ggplot(mapping = aes(x = `Th2 cells`, fill = Sex)) +
  geom_density(alpha = .5) +
  labs(title = "Th2 cell enrichment distributions", x = "Th2 cells xCell enrichment score", y = "Density", fill = "Gender") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 16))

# -------------------------------------------------------------------------------------------------------------------
## Saving ggplots

ggsave("Distri_R_S_CD8+_Tcm.png", plot = Distri_R_S_CD8plus_Tcm)
ggsave("Distri_R_S_CD8+_Tcells.png", plot = Distri_R_S_CD8plus_Tcells)
ggsave("Distri_R_S_Th1.png", plot = Distri_R_S_Th1_cells)
ggsave("Distri_R_S_Th2.png", plot = Distri_R_S_Th2_cells)
ggsave("Distri_G_CD8+_Tcm.png", plot = Distri_G_CD8plus_Tcm)
ggsave("Distri_G_CD8+_Tcells.png", plot = Distri_G_CD8plus_Tcells)
ggsave("Distri_G_Th1.png", plot = Distri_G_Th1_cells)
ggsave("Distri_G_Th2.png", plot = Distri_G_Th2_cells)
ggsave("Distri_R_S_gender_CD8+_Tcm.png", plot = Distri_R_S_gender_CD8plus_Tcm)
ggsave("Distri_R_S_gender_CD8+_Tcells.png", plot = Distri_R_S_gender_CD8plus_Tcells)
ggsave("Distri_R_S_gender_Th1.png", plot = Distri_R_S_gender_Th1_cells)
ggsave("Distri_R_S_gender_Th2.png", plot = Distri_R_S_gender_Th2_cells)

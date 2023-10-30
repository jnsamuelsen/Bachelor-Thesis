# Clearing global environment
rm(list = ls())

# -------------------------------------------------------------------------------------------------------------------
# Loading libraries

library(tidyverse)
library(ggplot2)
library(dplyr)
library(xCell)
library(matrixTests)
library(broom)
library(purrr)
library(modelr)
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
## 1. RRMS vs Healthy controls statistics

xCell_tibble_adj_R <- xCell_tibble_adj %>% 
  filter(Status == "R") %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)


xCell_tibble_adj_S <- xCell_tibble_adj %>% 
  filter(Status == "S") %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)


R_vs_S_wilcox <- col_wilcoxon_twosample(xCell_tibble_adj_R, xCell_tibble_adj_S)

## Boxplots

# CD8+ T-cells
RRMS_vs_HC_CD8plus_T_cells <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `CD8+ T-cells`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "CD8+ T-cell enrichment", x = "Health status", y = "CD8+ T-cells (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# CD8+ Tcm
RRMS_vs_HC_CD8plus_Tcm <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `CD8+ Tcm`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "CD8+ Tcm cell enrichment", x = "Health status", y = "CD8+ Tcm cells (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Th1 cells
RRMS_vs_HC_Th1_cells <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Th1 cells`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Th1 cell enrichment", x = "Health status", y = "Th1 cells (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Th2 cells
RRMS_vs_HC_Th2_cells <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Th2 cells`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Th2 cell enrichment", x = "Health status", y = "Th2 cells (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Smooth muscle
RRMS_vs_HC_Smooth_muscle_cells <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Smooth muscle`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Smooth muscle cell enrichment", x = "Health status", y = "Smooth muscle cells (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Neurons
RRMS_vs_HC_Neurons <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Neurons`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Neuron enrichment", x = "Health status", y = "Neurons (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Pericytes
RRMS_vs_HC_Pericytes <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Pericytes`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Pericyte enrichment", x = "Health status", y = "Pericytes (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Osteoblast
RRMS_vs_HC_Osteoblast <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Osteoblast`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Osteoblast enrichment", x = "Health status", y = "Osteoblasts (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# CMP
RRMS_vs_HC_CMP <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `CMP`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "CMP enrichment", x = "Health status", y = "CMP cells (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Platelets
RRMS_vs_HC_Platelets <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Platelets`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Platelet enrichment", x = "Health status", y = "Platelets (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Megakaryocytes
RRMS_vs_HC_Megakaryocytes <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Megakaryocytes`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Megakaryocyte enrichment", x = "Health status", y = "Megakaryocytes (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Neutrophils
RRMS_vs_HC_Neutrophils <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Neutrophils`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Neutrophil enrichment", x = "Health status", y = "Neutrophils (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# CD8+ Tem
RRMS_vs_HC_CD8plus_Tem <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `CD8+ Tem`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "CD8+ Tem cell enrichment", x = "Health status", y = "CD8+ Tem cells (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Erythrocytes
RRMS_vs_HC_Erythrocytes <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Erythrocytes`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Erythrocyte enrichment", x = "Health status", y = "Erythrocytes (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))

# Fibroblasts
RRMS_vs_HC_Fibroblasts <- xCell_tibble_adj %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Fibroblasts`,
                       fill = Status)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Status), size = 6, label.x.npc = "left") +
  labs(title = "Fibroblast enrichment", x = "Health status", y = "Fibroblasts (xCell enrichment)", fill = "Health status") +
  scale_x_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  scale_fill_discrete(labels = c("R" = "RRMS", "S" = "HC")) +
  theme(text = element_text(size = 20))


# -------------------------------------------------------------------------------------------------------------------
## 2. Active vs inactive MS patients on fingolimod statistics

# Grouping data by stability
xCell_tibble_adj_G_unstabile <- xCell_tibble_adj %>% 
  filter(Status == "G", Stability == 1) %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

xCell_tibble_adj_G_stabile <- xCell_tibble_adj %>% 
  filter(Status == "G", Stability == 0) %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

stabile_vs_unstabile_wilcox<- col_wilcoxon_twosample(xCell_tibble_adj_G_stabile, xCell_tibble_adj_G_unstabile)

## Boxplots

# Renaming Stability column values to 0 ~ stabile, 1 ~ unstabile
xCell_tibble_adj <- xCell_tibble_adj %>% 
  mutate(Stability = case_when(Stability == 0 ~ "Stabile",
                               Stability == 1 ~ "Unstabile"))
# Smooth muscle
G_stabile_vs_unstabile_Smooth_muscle_cells <- xCell_tibble_adj %>% 
  filter(Status == "G") %>% 
  ggplot(mapping = aes(x = Stability, y = `Smooth muscle`,
                       fill = Stability)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Stability), size = 6, label.x.npc = "center") +
  labs(title = "Smooth muscle cell enrichment", x = "MS patient stability", y = "Smooth muscle cells (xCell enrichment)", fill = "MS patient stability") +
  scale_fill_brewer(palette = "Dark2") +
  theme(text = element_text(size = 20))

# pDC
G_stabile_vs_unstabile_pDC_cells <- xCell_tibble_adj %>% 
  filter(Status == "G") %>% 
  ggplot(mapping = aes(x = Stability, y = `pDC`,
                       fill = Stability)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Stability), size = 6, label.x.npc = "center") +
  labs(title = "pDC cell enrichment", x = "MS patient stability", y = "pDC cells (xCell enrichment)", fill = "MS patient stability") +
  scale_fill_brewer(palette = "Dark2") +
  theme(text = element_text(size = 20))

# Class_switched_memory_Bcells
G_stabile_vs_unstabile_class_switched_memory_Bcells <- xCell_tibble_adj %>% 
  filter(Status == "G") %>% 
  ggplot(mapping = aes(x = Stability, y = `Class-switched memory B-cells`,
                       fill = Stability)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Stability), size = 6, label.x.npc = "center") +
  labs(title = "Class-switched memory B-cell enrichment", x = "MS patient stability", y = "Class-switched memory B-cells (xCell enrichment)", fill = "MS patient stability") +
  scale_fill_brewer(palette = "Dark2") +
  theme(text = element_text(size = 20),
        plot.title = element_text(size = 18),
        axis.title.y = element_text(size = 17.5)) 



# -------------------------------------------------------------------------------------------------------------------
## 3. RRMS males vs RRMS females

xCell_tibble_age_adj_males <- xCell_tibble_age_adj %>% 
  filter(Status == "R", Sex == 1) %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

xCell_tibble_age_adj_females <- xCell_tibble_age_adj %>% 
  filter(Status == "R", Sex == 0) %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

R_males_vs_females_wilcox <- col_wilcoxon_twosample(xCell_tibble_age_adj_males, xCell_tibble_age_adj_females)

## Boxplots

# Renaming Sex column values to 0 ~ Female, 1 ~ Male
xCell_tibble <- xCell_tibble %>% 
  mutate(Sex = case_when(Sex == 0 ~ "Female",
                               Sex == 1 ~ "Male"))

# CMP
RRMS_males_vs_females_CMP <- xCell_tibble %>% 
  filter(Status == "R") %>% 
  ggplot(mapping = aes(x = Sex, y = `CD8+ Tem`,
                       fill = Sex)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Sex), size = 6, label.x.npc = "center") +
  labs(title = "CMP enrichment", x = "Gender", y = "CMP cells (xCell enrichment)", fill = "Gender") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 20))

# Melanocytes
RRMS_males_vs_females_Melanocytes <- xCell_tibble %>% 
  filter(Status == "R") %>% 
  ggplot(mapping = aes(x = Sex, y = `Melanocytes`,
                       fill = Sex)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Sex), size = 6, label.x.npc = "center") +
  labs(title = "Melanocyte enrichment", x = "Gender", y = "Melanocytes (xCell enrichment)", fill = "Gender") +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 20))

# CMP
RRMS_males_vs_females_CMP <- xCell_tibble_age_adj %>% 
  filter(Status == "R") %>% 
  ggplot(mapping = aes(x = Sex, y = `CMP`,
                       fill = Sex)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Sex), size = 6, label.x.npc = "left") +
  labs(title = "CMP enrichment", x = "Gender", y = "CMP cells (xCell enrichment)", fill = "Gender", size = 20) +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 20))

# CD8+ Tem
RRMS_males_vs_females_CD8plus_Tem <- xCell_tibble_age_adj %>% 
  filter(Status == "R") %>% 
  ggplot(mapping = aes(x = Sex, y = `CD8+ Tem`,
                       fill = Sex)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Sex), size = 6, label.x.npc = "left") +
  labs(title = "CD8+ Tem cell enrichment", x = "Gender", y = "CD8+ Tem cells (xCell enrichment)", fill = "Gender", size = 20) +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 20))

# NK cells
RRMS_males_vs_females_NK_cells <- xCell_tibble_age_adj %>% 
  filter(Status == "R") %>% 
  ggplot(mapping = aes(x = Sex, y = `NK cells`,
                       fill = Sex)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Sex), size = 6, label.x.npc = "left") +
  labs(title = "NK cell enrichment", x = "Gender", y = "NK cells (xCell enrichment)", fill = "Gender", size = 20) +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 20))

# Boxplots for CD8+ Tem and NK cells doesn't show anything
# Calculating some summary statistics for CD8+ Tem and NK cells

# CD8+ Tem
RRMS_males_vs_females_CD8plus_Tem_table <- xCell_tibble_age_adj %>% 
  filter(Status == "R") %>% 
  group_by(Sex) %>%
  summarise(
    count = n(),
    mean = mean(`CD8+ Tem`, na.rm = TRUE),
    median = median(`CD8+ Tem`, na.rm = TRUE),
    IQR = IQR(`CD8+ Tem`, na.rm = TRUE))

# NK cells
RRMS_males_vs_females_NK_cells_table <-xCell_tibble_age_adj %>% 
  filter(Status == "R") %>% 
  group_by(Sex) %>%
  summarise(
    count = n(),
    mean = mean(`NK cells`, na.rm = TRUE),
    median = median(`NK cells`, na.rm = TRUE),
    IQR = IQR(`NK cells`, na.rm = TRUE))

# -------------------------------------------------------------------------------------------------------------------
## Saving ggplots

ggsave("R_S_CD8+_Tcells.png", plot = RRMS_vs_HC_CD8plus_T_cells)
ggsave("R_S_CD8+_Tcm.png", plot = RRMS_vs_HC_CD8plus_Tcm)
ggsave("R_S_Th1.png", plot = RRMS_vs_HC_Th1_cells)
ggsave("R_S_Th2.png", plot = RRMS_vs_HC_Th2_cells)
ggsave("R_S_Smooth_muscle.png", plot = RRMS_vs_HC_Smooth_muscle_cells)
ggsave("R_S_Neurons.png", plot = RRMS_vs_HC_Neurons)
ggsave("R_S_Pericytes.png", plot = RRMS_vs_HC_Pericytes)
ggsave("R_S_Osteoblast.png", plot = RRMS_vs_HC_Osteoblast)
ggsave("R_S_CMP.png", plot = RRMS_vs_HC_CMP)
ggsave("R_S_Platelets.png", plot = RRMS_vs_HC_Platelets)
ggsave("R_S_Megakaryocytes.png", plot = RRMS_vs_HC_Megakaryocytes)
ggsave("R_S_Neutrophils.png", plot = RRMS_vs_HC_Neutrophils)
ggsave("R_S_CD8+_Tem.png", plot = RRMS_vs_HC_CD8plus_Tem)
ggsave("R_S_Erythrocytes.png", plot = RRMS_vs_HC_Erythrocytes)
ggsave("R_S_Fibroblasts.png", plot = RRMS_vs_HC_Fibroblasts)
ggsave("G_S_U_Smooth_muscle.png", plot = G_stabile_vs_unstabile_Smooth_muscle_cells)
ggsave("G_S_U_pDC.png", plot = G_stabile_vs_unstabile_pDC_cells)
ggsave("G_S_U_Class_switched_memory_Bcells.png", plot = G_stabile_vs_unstabile_class_switched_memory_Bcells)
ggsave("R_M_F_CMP.png", plot = RRMS_males_vs_females_CMP)
ggsave("R_M_F_CD8+_Tem.png", plot = RRMS_males_vs_females_CD8plus_Tem)
ggsave("R_M_F_NK_cells.png", plot = RRMS_males_vs_females_NK_cells)

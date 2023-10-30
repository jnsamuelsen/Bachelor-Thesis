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

# Loading data into R
load("C:/Users/jeppe/Desktop/Bachelorprojekt/Data/MS data/data.RData")

# Running xCell on expression matrix
xCell_matrix <- xCellAnalysis(matrix_exprs)

# Transposing matrix
xCell_matrix <- t(xCell_matrix)

# Creating a tibble for the data
xCell_tibble <- xCell_matrix %>% 
  as_tibble()

# Adding the filename column to the tibble
xCell_tibble <- xCell_tibble %>%
  mutate(Filename = rownames(xCell_matrix))

# Joining the pheno data to the xCell tibble
xCell_tibble <- full_join(xCell_tibble, pheno, by = "Filename")

# Renaming columns
xCell_tibble <- xCell_tibble %>% 
  rename(Status = R.RRMS..S.HC..G.Gilenya.behandlet.,
         Sex = sex.male1,
         Past_treatment = Tidligere.Beh,
         Stability = X0Stabil_1ustabil)

# -------------------------------------------------------------------------------------------------------------------
# 1. RRMS vs Healthy controls
# Disse dele (både 1. og 2. kan omskrives til bare at bruger "filter" og "select", og så lave tests direkte evt.)
# Mann-Whitney U tests (Wilcoxon tests)

# Creating separate tibbles for enrichment data of R (RRMS) and S (HC) in the "Treatment" column

xCell_tibble_R <- xCell_tibble %>% 
  filter(Status == "R") %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)


xCell_tibble_S <- xCell_tibble %>% 
  filter(Status == "S") %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

# Performing the statistical test of independence
R_vs_S_wilcox <- col_wilcoxon_twosample(xCell_tibble_R, xCell_tibble_S)

# Now we do the same test for only males
xCell_tibble_males_R <- xCell_tibble %>% 
  filter(Sex == 1, Status == "R") %>%
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

xCell_tibble_males_S <- xCell_tibble %>% 
  filter(Sex == 1, Status == "S") %>%
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

# Performing the statistical test of independence
males_R_vs_S_wilcox <- col_wilcoxon_twosample(xCell_tibble_males_R, xCell_tibble_males_S)


# Now we do the same test for only females
xCell_tibble_females_R <- xCell_tibble %>% 
  filter(Sex == 0, Status == "R") %>%
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

xCell_tibble_females_S <- xCell_tibble %>% 
  filter(Sex == 0, Status == "S") %>%
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

# Performing the statistical test of independence
females_R_vs_S_wilcox <- col_wilcoxon_twosample(xCell_tibble_females_R, xCell_tibble_females_S)

# Grouping data by treatment and sex, then nest the data
by_gender <- xCell_tibble %>% 
  group_by(Status, Sex) %>% 
  nest()

# Defining Fibroblasts model
Fibroblasts_model <- function(df) {
  lm(Fibroblasts ~ Age, data = df)
}

# Applying models to the data
models <- map(by_gender$data, Fibroblasts_model)

# Adding the model column to the by_gender tibble
by_gender <- by_gender %>% 
  mutate(model = map(data, Fibroblasts_model))

# Adding residuals to the by_gender tibble
by_gender <- by_gender %>% 
  mutate(
    resids = map2(data, model, add_residuals)
  )

# Unnesting resids for plotting
resids <- unnest(by_gender, resids)

resids %>% 
  ggplot(mapping = aes(Age, resid, group = Status)) +
  geom_line() +
  facet_wrap(~Sex)


fit <- lm(`CD8+ T-cells` ~ Age, data = xCell_tibble)

summary(fit) %>% map_dfr(glance)
summary(fit) %>% map_dfr(tidy)

models <- xCell_tibble %>% 
  pivot_longer(
    cols = aDC:Tregs,
    names_to = "y_name",
    values_to = "y_value"
  ) %>%
  split(.$y_name) %>%
  map(~lm(y_value ~ Age, data = .)) %>%
  tibble(
    dvsub = names(.),
    untidied = .
  ) %>%
  mutate(tidy = map(untidied, broom::tidy)) %>%
  unnest(tidy) 

### Boxplots showing Mann-Whitney U test p-values (facetted on gender/sex) ###

# Creating name label vector for facet variable (Sex)
Sex_labs <- c("Females", "Males")
names(Sex_labs) <- c(0, 1)

# CD8+ T-cells
xCell_tibble %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `CD8+ T-cells`,
                                            fill = Status)) +
  geom_boxplot() +
  facet_wrap(~ Sex,
             labeller = labeller(Sex = Sex_labs)) +
  stat_compare_means(aes(group = Status)) +
  labs(title = "CD8+ T-cell enrichment", x = "Patient status", y = "CD8+ T-cells (xCell enrichment)", fill = "Patient status")

# Th1 cells
xCell_tibble %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Th1 cells`,
                       fill = Status)) +
  geom_boxplot() +
  facet_wrap(~ Sex,
             labeller = labeller(Sex = Sex_labs)) +
  stat_compare_means(aes(group = Status)) +
  labs(title = "Th1 cell enrichment", x = "Patient status", y = "Th1 cells (xCell enrichment)", fill = "Patient status")

# CD8+ Tcm
xCell_tibble %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `CD8+ Tcm`,
                       fill = Status)) +
  geom_boxplot() +
  facet_wrap(~ Sex,
             labeller = labeller(Sex = Sex_labs)) +
  stat_compare_means(aes(group = Status)) +
  labs(title = "CD8+ Tcm enrichment", x = "Patient status", y = "CD8+ Tcm cells (xCell enrichment)", fill = "Patient status")

# Th2 cells
xCell_tibble %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Th2 cells`,
                       fill = Status)) +
  geom_boxplot() +
  facet_wrap(~ Sex,
             labeller = labeller(Sex = Sex_labs)) +
  stat_compare_means(aes(group = Status)) +
  labs(title = "Th2 cell enrichment", x = "Patient status", y = "Th2 cells (xCell enrichment)", fill = "Patient status")

# Tregs
xCell_tibble %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Tregs`,
                       fill = Status)) +
  geom_boxplot() +
  facet_wrap(~ Sex,
             labeller = labeller(Sex = Sex_labs)) +
  stat_compare_means(aes(group = Status)) +
  labs(title = "Treg cell enrichment", x = "Patient status", y = "Treg cells (xCell enrichment)", fill = "Patient status")

# CMP
xCell_tibble %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `CMP`,
                       fill = Status)) +
  geom_boxplot() +
  facet_wrap(~ Sex,
             labeller = labeller(Sex = Sex_labs)) +
  stat_compare_means(aes(group = Status)) +
  labs(title = "CMP cell enrichment", x = "Patient status", y = "CMP cells (xCell enrichment)", fill = "Patient status")

# Plasma cells
xCell_tibble %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Plasma cells`,
                       fill = Status)) +
  geom_boxplot() +
  facet_wrap(~ Sex,
             labeller = labeller(Sex = Sex_labs)) +
  stat_compare_means(aes(group = Status)) +
  labs(title = "Plasma cell enrichment", x = "Patient status", y = "Plasma cells (xCell enrichment)", fill = "Patient status")

# Megakaryocytes
xCell_tibble %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Megakaryocytes`,
                       fill = Status)) +
  geom_boxplot() +
  facet_wrap(~ Sex,
             labeller = labeller(Sex = Sex_labs)) +
  stat_compare_means(aes(group = Status)) +
  labs(title = "Megakaryocyte enrichment", x = "Patient status", y = "Megakaryocyte cells (xCell enrichment)", fill = "Patient status")

# Sebocytes
xCell_tibble %>% 
  filter(Status == "R" | Status == "S") %>% 
  ggplot(mapping = aes(x = Status, y = `Sebocytes`,
                       fill = Status)) +
  geom_boxplot() +
  facet_wrap(~ Sex,
             labeller = labeller(Sex = Sex_labs)) +
  stat_compare_means(aes(group = Status)) +
  labs(title = "Sebocyte enrichment", x = "Patient status", y = "Sebocyte cells (xCell enrichment)", fill = "Patient status")


# -------------------------------------------------------------------------------------------------------------------
# 2. Active vs inactive MS patients on fingolimod Status

# Grouping data by stability
xCell_tibble_Active_Fingolimod <- xCell_tibble %>% 
  filter(Status == "G", Stability == 1) %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

xCell_tibble_Inactive_Fingolimod <- xCell_tibble %>% 
  filter(Status == "G", Stability == 0) %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

# Performing the statistical test of independence for active vs inactive MS patients on fingolimod
Fingolimod_Active_vs_Inactive_wilcox <- col_wilcoxon_twosample(xCell_tibble_Active_Fingolimod, xCell_tibble_Inactive_Fingolimod)

# Now do the same test for males
xCell_tibble_males_Active_Fingolimod <- xCell_tibble %>% 
  filter(Status == "G", Stability == 1, Sex == 1) %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

xCell_tibble_males_Inactive_Fingolimod <- xCell_tibble %>% 
  filter(Status == "G", Stability == 0, Sex == 1) %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

# Performing the statistical test of independence for active vs inactive male MS patients on fingolimod
Fingolimod_males_Active_vs_Inactive_wilcox <- col_wilcoxon_twosample(xCell_tibble_males_Active_Fingolimod, xCell_tibble_males_Inactive_Fingolimod)

# Now do the same test for females
xCell_tibble_females_Active_Fingolimod <- xCell_tibble %>% 
  filter(Status == "G", Stability == 1, Sex == 0) %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

xCell_tibble_females_Inactive_Fingolimod <- xCell_tibble %>% 
  filter(Status == "G", Stability == 0, Sex == 0) %>% 
  select(-ImmuneScore, -StromaScore, -MicroenvironmentScore, -Filename, -Status, -Sex, -Age, -Affy.Run_nr, -Stability, -Past_treatment)

# Performing the statistical test of independence for active vs inactive female MS patients on fingolimod
Fingolimod_females_Active_vs_Inactive_wilcox <- col_wilcoxon_twosample(xCell_tibble_females_Active_Fingolimod, xCell_tibble_females_Inactive_Fingolimod)

# Renaming stability column values to 0 ~ stabile, 1 ~ unstabile
xCell_tibble <- xCell_tibble %>% 
  mutate(Stability = case_when(Stability == 0 ~ "Stabile",
                               Stability == 1 ~ "Unstabile"))

# Boxplots for active vs inactive patients facetted on earlier Status
# CD8+ T-cells
xCell_tibble %>% 
  filter(Status == "G") %>% 
  ggplot(mapping = aes(x = Stability, y = `CD8+ T-cells`,
                       fill = Stability)) +
  geom_boxplot() +
  facet_wrap(~ Past_treatment) +
  stat_compare_means(aes(group = Stability)) +
  labs(title = "CD8+ T-cell enrichment", x = "Patient MS stability", y = "CD8+ T-cells (xCell enrichment)", fill = "Patient MS stability")

# Smooth muscle cells
xCell_tibble %>% 
  filter(Status == "G") %>% 
  ggplot(mapping = aes(x = Stability, y = `Smooth muscle`,
                       fill = Stability)) +
  geom_boxplot() +
  facet_wrap(~ Past_treatment) +
  stat_compare_means(aes(group = Stability)) +
  labs(title = "Smooth muscle cells enrichment", x = "Patient MS stability", y = "Smooth muscle cells (xCell enrichment)", fill = "Patient MS stability")

# Performing the statistical test of independence for RRMS male vs RRMS female patients
RRMS_females_vs_RRMS_males_wilcox <- col_wilcoxon_twosample(xCell_tibble_females_R, xCell_tibble_males_R)


## Save tables ##
print(xtable(R_vs_S_wilcox), type = "latex", file = "R_vs_S_stats.tex")

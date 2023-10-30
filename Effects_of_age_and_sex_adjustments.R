# Clearing global environment
rm(list = ls())

# -------------------------------------------------------------------------------------------------------------------
# Loading libraries

library(tidyverse)
library(ggplot2)
library(dplyr)
library(broom)
library(purrr)
library(modelr)
library(tidyr)
library(gridExtra)

# Loading data into R
load("C:/Users/jeppe/Desktop/Bachelorprojekt/Data/MS data/data_updated.RData")

pheno <- pheno %>% 
  mutate(sex.male1 = case_when(sex.male1 == 0 ~ "Female",
                         sex.male1 == 1 ~ "Male"))

## Find most severe example of sex adjustment
r2s <- c()
cor <- c()
slopes <- c()
for(j in 1:nrow(matrix_exprs)) {
  fit <- lm(matrix_exprs_sex_adj[j,] ~ matrix_exprs[j,])
  r2s <- append(r2s, summary(fit)$adj.r.squared)
  cor <- append(cor, cor(matrix_exprs_sex_adj[j,], matrix_exprs[j,]))
  slopes <- append(slopes, fit$coef[[2]])
}
rownames(matrix_exprs)[which.min(r2s)]
rownames(matrix_exprs)[which.min(cor)]
rownames(matrix_exprs)[which.max(1-abs(slopes))]
gene <- which.min(r2s)


## Plot most severe case
fit <- lm(matrix_exprs[gene,] ~ pheno$sex.male1)
data <- data.frame(sex=as.factor(pheno$sex.male1), expr=matrix_exprs[gene,])
g1_gender <- ggplot(data, aes(x=sex, y=expr, fill=sex)) +
  geom_boxplot() +
  annotate(geom="text", x=2, y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=0.5, size = 5) +
  labs(title = paste0("Gender vs expression of ", rownames(matrix_exprs)[gene]), fill = "Gender", x = "Gender", y = "Expression") +
  theme(text = element_text(size = 20))

fit <- lm(matrix_exprs_sex_adj[gene,] ~ matrix_exprs[gene,])
data <- data.frame(expr_adj=matrix_exprs_sex_age_adj[gene,], expr=matrix_exprs[gene,], sex=as.factor(pheno$sex.male1))
g2_gender <- ggplot(data, aes(x=expr_adj, y=expr)) +
  geom_point(aes(color=sex)) +
  geom_smooth(method=lm) +
  annotate(geom="text", x=max(data$expr_adj), y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=2, size = 5) +
  labs(title = paste0("Gender adjusted expression vs expression of ", rownames(matrix_exprs_sex_age_adj)[gene]), color = "Gender", x = "Gender adjusted expression", y = "Expression") +
  theme(text = element_text(size = 19))

fit <- lm(matrix_exprs_sex_adj[gene,] ~ pheno$sex.male1)
data <- data.frame(sex=as.factor(pheno$sex.male1), expr=matrix_exprs_sex_adj[gene,])
g3_gender <- ggplot(data, aes(x=sex, y=expr, fill=sex)) +
  geom_boxplot() +
  annotate(geom="text", x=2, y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=1.5, size = 5) +
  labs(title = paste0("Gender vs gender adjusted expression of ", rownames(matrix_exprs)[gene]), color = "Gender", x = "Gender", y = "Gender adjusted expression", fill = "Gender") +
  theme(text = element_text(size = 20))

# grid.arrange(g1,g2,g3, ncol=3, nrow=1)


## Find most severe case of age adjustment
r2s <- c()
cor <- c()
slopes <- c()
for(j in 1:nrow(matrix_exprs)) {
  fit <- lm(matrix_exprs_sex_age_adj[j,] ~ matrix_exprs_sex_adj[j,])
  r2s <- append(r2s, summary(fit)$adj.r.squared)
  cor <- append(cor, cor(matrix_exprs_sex_age_adj[j,], matrix_exprs_sex_adj[j,]))
  slopes <- append(slopes, fit$coef[[2]])
}
rownames(matrix_exprs)[which.min(r2s)]
rownames(matrix_exprs)[which.min(cor)]
rownames(matrix_exprs)[which.max(1-abs(slopes))]
gene <- which.min(r2s)


## Plot most severe case of age adjustment
fit <- lm(matrix_exprs_sex_adj[gene,] ~ pheno$Age)
data <- data.frame(age=pheno$Age, expr=matrix_exprs_sex_adj[gene,])
g1_age <- ggplot(data, aes(x=age, y=expr)) +
  geom_point() +
  geom_smooth(method=lm) +
  annotate(geom="text", x=max(data$age), y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=1, size = 5) +
  ggtitle(paste0("Age vs expression of ", rownames(matrix_exprs)[gene])) +
  labs(x = "Age", y = "Expression") +
  theme(text = element_text(size = 20))

fit <- lm(matrix_exprs_sex_age_adj[gene,] ~ matrix_exprs_sex_adj[gene,])
data <- data.frame(expr_adj=matrix_exprs_sex_age_adj[gene,], expr=matrix_exprs_sex_adj[gene,])
g2_age <- ggplot(data, aes(x=expr_adj, y=expr)) +
  geom_point() +
  geom_smooth(method=lm) +
  annotate(geom="text", x=max(data$expr_adj), y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=2, size = 5) +
  ggtitle(paste0("Age adjusted expression vs expression of ", rownames(matrix_exprs_sex_age_adj)[gene])) +
  labs(x = "Age adjusted expression", y = "Expression") +
  theme(text = element_text(size = 19))

fit <- lm(matrix_exprs_sex_age_adj[gene,] ~ pheno$Age)
data <- data.frame(age=pheno$Age, expr=matrix_exprs_sex_age_adj[gene,])
g3_age <- ggplot(data, aes(x=age, y=expr)) +
  geom_point() +
  geom_smooth(method=lm) +
  annotate(geom="text", x=max(data$age), y=max(data$expr), label = paste0("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "\nSlope = ",signif(fit$coef[[2]], 5)), hjust=1, size = 5) +
  ggtitle(paste0("Age vs age adjusted expression of ", rownames(matrix_exprs)[gene])) +
  labs(x = "Age", y = "Age adjusted expression") +
  theme(text = element_text(size = 20))

# grid.arrange(g1,g2,g3, ncol=3, nrow=1)

# -------------------------------------------------------------------------------------------------------------------
## Saving ggplots
ggsave("XIST_Gender_vs_Expr.png", plot = g1_gender)
ggsave("XIST_Gender_adj_expr_vs_Expr.png", plot = g2_gender)
ggsave("XIST_Gender_vs_Gender_adj_expr.png", plot = g3_gender)
ggsave("NRCAM_Age_vs_Expr.png", plot = g1_age)
ggsave("NRCAM_Age_adj_expr_vs_Expr.png", plot = g2_age)
ggsave("NRCAM_Age_vs_age_adj_expr.png", plot = g3_age)

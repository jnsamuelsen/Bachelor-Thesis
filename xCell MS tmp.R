library(tidyverse)
library(ggplot2)
library(dplyr)
library(xCell)
library(matrixTests)
library(broom)

ggplot(data = xCell_tibble, mapping = aes(x = "", y = Neutrophils, colour = sex.male1)) +
  geom_boxplot() +
  facet_wrap(~sex.male1)

df1 <- data.frame(
  g = runif(50), 
  pair = sample(x = c("A", "B", "C"), size = 50, replace = TRUE), 
  V1 = runif(50), 
  V2 = runif(50), 
  V3 = runif(50), 
  V4 = runif(50), 
  V5 = runif(50),
  stringsAsFactors = FALSE
)

df2 <- df1 %>% 
  as_tibble %>% 
  gather(key = "column", value = "value", V1:V5) %>%       # first set the data in long format
  nest(g, value) %>%                                       # now nest the dependent and independent factors
  mutate(model = map(data, ~lm(g ~ value, data = .))) %>%  # fit the model using purrr
  mutate(tidy_model = map(model, tidy)) %>%                # clean the model output with broom
  select(-data, -model) %>%                                # remove the "untidy" parts
  unnest()  



# Linear models to test for age / cell enrichment correlations
age_models_tibble <- xCell_tibble %>% 
  gather(key = "column", value = "value", aDC:Tregs) %>% 
  nest(Age, value) %>% 
  mutate(model = map(data, ~lm(value ~ Age, data = .))) %>% 
  mutate(tidy_model = map(model, tidy)) %>%
  select(-data, -model) %>%   
  unnest()




# -------------------------------------------------------------------------------------------------------------------
# 2. Active vs inactive MS patients on fingolimod treatment

# Grouping data by stability
by_stability <- xCell_tibble %>% 
  filter(Treatment == "G") %>% 
  pivot_longer(
    cols = aDC:Tregs,
    names_to = "cell_type",
    values_to = "enrichment_value"
  ) %>%
  group_by(Stability, cell_type)

# Summarizing mean enrichment values
summarise(by_stability, cell_enrichment = mean(enrichment_value, na.rm = TRUE))

by_stability %>% 
  mutate(enrichment_value = as.numeric(enrichment_value)) %>% 
  drop_na() %>% 
  group_by(cell_type) %>% 
  wilcox.test(enrichment_value ~ Stability)
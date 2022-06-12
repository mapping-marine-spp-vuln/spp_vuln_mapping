library(tidyverse)

n_spp <- 1000
n_cells <- 1000

df <- crossing(spp = 1:n_spp, cell = 1:n_cells) %>%
  mutate(impact = rnorm(.3, .1, n = n()) %>% abs(),
         effect_convex = impact^(2),
         effect_concave = impact^(0.5))

ggplot(df %>% sample_n(1000), 
       aes(x = impact, y = effect_convex)) + 
  geom_point()
ggplot(df %>% sample_n(1000), 
       aes(x = impact, y = effect_concave)) + 
  geom_point()

sum_df <- df %>%
  group_by(spp) %>%
  summarize(mean_eff_conv = mean(effect_convex),
            mean_eff_conc = mean(effect_concave),
            mean_imp = mean(impact),
            sd_imp   = sd(impact),
            cov_imp  = sd(impact) / mean(impact))

lm(mean_eff_conv ~ mean_imp + sd_imp, data = sum_df)
lm(mean_eff_conc ~ mean_imp + sd_imp, data = sum_df)
lm(mean_eff_conv ~ mean_imp + cov_imp, data = sum_df)
lm(mean_eff_conc ~ mean_imp + cov_imp, data = sum_df)

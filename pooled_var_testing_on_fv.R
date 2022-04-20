### does pooled var work for functional vulnerability/FEs?

df <- data.frame(
  fe = sample(1:20, size = 100, replace = TRUE),
  v = rnorm(mean = .25, sd = .25, n = 100)) %>%
  mutate(v = case_when(v < 0 ~ 0,
                       v > 1 ~ 1,
                       TRUE ~ v)) %>%
  group_by(fe) %>%
  mutate(n_fe = n(),
         fv = .5^(n_fe - 1)) %>%
  ungroup()
fv_mean_df <- df %>%
  group_by(fe) %>%
  summarize(mean_v = mean(v),
            sdev_v = sd(v),
            fv = first(fv),
            n_fe = first(n_fe),
            .groups = 'drop')

fv_sdev_overall <- fv_mean_df %>%
  summarize(mean_v_all = Hmisc::wtd.mean(mean_v, weights = fv), 
            sd_v_all   = Hmisc::wtd.var(mean_v, weights = fv))
  
fv_a <- fv_mean_df %>% 
  filter(fe < 6) %>%
  summarize(mean_v = mean(mean_v), 
            sd_v = sd(mean_v),
            sum_fv = sum(fv))
fv_b <- fv_mean_df %>% 
  filter(between(fe, 6, 15)) %>%
  summarize(mean_v = mean(mean_v), 
            sd_v = sd(mean_v),
            sum_fv = sum(fv))
fv_c <- fv_mean_df %>%
  filter(fe > 15) %>%
  summarize(mean_v = mean(mean_v), 
            sd_v = sd(mean_v),
            sum_fv = sum(fv))

fv_chunks <- bind_rows(fv_a, fv_b, fv_c)


pooled_var <- function(x_bar, y_bar, s_x, s_y, n_x, n_y) {
  ### convert std dev to var
  var_x <- ifelse(is.na(s_x), 0, s_x^2)
  var_y <- ifelse(is.na(s_y), 0, s_y^2)
  
  var_xy_clean <- ((n_x - 1)*var_x + (n_y - 1)*var_y) / (n_x + n_y - 1)
  var_xy_error <- (n_x * n_y) * (x_bar - y_bar)^2 / ((n_x + n_y)*(n_x + n_y - 1))
  
  return(var_xy_clean + var_xy_error)
}

iterated_pooled_var <- function(mean_vec, sdev_vec, n_vec, flag = FALSE) {
  if(!all.equal(length(mean_vec), length(sdev_vec), length(n_vec))) {
    stop('Mean, std dev, and n vectors must all be equal length!')
  }
  if(length(mean_vec) == 1) {
    warning('Only one element - no need for pooled variance!')
    return(sdev_vec[1]^2)
  }
  ### initialize values for first in list
  mean_x <- mean_vec[1]; s_x <- sdev_vec[1]; n_x <- n_vec[1]
  for(i in 2:length(mean_vec)) { ## i <- 2
    if(flag) message('  ... processing iteration ', i - 1, '...')
    
    mean_y <- mean_vec[i]
    s_y    <- sdev_vec[i]
    n_y    <- n_vec[i]
    var_out <- pooled_var(x_bar = mean_x, y_bar = mean_y, 
                          n_x = n_x, n_y = n_y, 
                          s_x = s_x, s_y = s_y)
    
    ### set up values for next iteration
    mean_x <- (mean_x * n_x + mean_y * n_y) / (n_x + n_y)
    s_x <- sqrt(var_out)
    n_x <- n_x + n_y
  }
  return(var_out)
}


pooled_var_results <- iterated_pooled_var(mean_vec = fv_chunks$mean_v,
                                          sdev_vec = fv_chunks$sd_v,
                                          n_vec = fv_chunks$sum_fv)

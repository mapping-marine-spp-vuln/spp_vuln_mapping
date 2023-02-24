#######################################################.
#### Helper functions for FishBase and SeaLifeBase ####
#######################################################.

### Get vulnerability traits by species
get_vuln_traits <- function() {
  x <- data.table::fread(here('_data/traits_vulnerability/spp_traits_valid.csv'))
  return(x)
}

### get traits from FishBase/SeaLifeBase, keeping or dropping columns as needed
get_fb_slb <- function(fxn = species, 
                       keep_cols = NULL, keep_fxn = contains, 
                       drop_cols = NULL, drop_fxn = all_of) {
  fb <- fxn(server = 'fishbase') %>%
    janitor::clean_names() %>%
    mutate(db = 'fb')
  slb <- fxn(server = 'sealifebase') %>%
    janitor::clean_names() %>%
    mutate(db = 'slb')
  if(!is.null(keep_cols)) {
    fb <- fb %>%
      dplyr::select(spec_code, species, db, keep_fxn(keep_cols))
    slb <- slb %>%
      dplyr::select(spec_code, species, db, keep_fxn(keep_cols))
  }
  if(!is.null(drop_cols)) {
    fb <- fb %>%
      dplyr::select(-drop_fxn(drop_cols))
    slb <- slb %>%
      dplyr::select(-drop_fxn(drop_cols))
  }
  ### sometimes class of a column from one server mismatches that of the other.
  ### This should only be a problem for character vs. numeric - i.e.,
  ### not a problem for numeric vs. integer
  cols_check <- data.frame(
    type_fb = sapply(fb, class),
    col = names(fb)) %>%
    inner_join(data.frame(
      type_slb = sapply(slb, class),
      col = names(slb))) %>%
    filter(type_fb != type_slb) %>%
    filter(type_fb == 'character' | type_slb == 'character')
  if(nrow(cols_check) > 0) {
    message('Conflicting column types; coercing all to character:')
    message(cols_check)
    fb <- fb %>%
      mutate(across(all_of(cols_check$col), as.character))
    slb <- slb %>%
      mutate(across(all_of(cols_check$col), as.character))
  }
  
  return(bind_rows(fb, slb) %>% mutate(species = tolower(species)))
}

###################################################.
#### function for downstream direct match fill ####
###################################################.

assign_repr_level <- function(df) {
  # message('Assigning representative ranks to species-level info...')
  # repres_fish <- readxl::read_excel(here('_raw/xlsx/fish_traits_add_bp.xlsx'), 
  #                                   sheet = 'repres_expanded')
  # write_csv(repres_fish, here('_raw/spp_gp_repres_fish.csv'))
  repres_fish <- read_csv(here('_raw/spp_gp_repres_fish.csv'), show_col_types = FALSE)
  # repres_df   <- readxl::read_excel(here('_raw/xlsx/spp_gp_representativeness.xlsx')) %>%
  #   janitor::clean_names() %>%
  #   select(taxon, species, repres = representative_rank)
  # write_csv(repres_df, here('_raw/spp_gp_repres_all_else.csv'))
  repres_df   <- read_csv(here('_raw/spp_gp_repres_all_else.csv'), show_col_types = FALSE) %>%
    bind_rows(repres_fish) %>%
    mutate(species = tolower(species),
           repres = tolower(repres)) %>%
    filter(!is.na(repres) & repres != 'species')
  
  if(any(c('phylum', 'class', 'order', 'family', 'genus') %in% names(df))) {
    stop('get rid of taxonomic ranks higher than species!')
  }
  spp_all_wide <- assemble_worms(aspect = 'wide')
  rank_names_long <- assemble_worms(aspect = 'long') %>%
    select(-spp) %>%
    distinct()
  
  ### create new dataframe including all representative species, with
  ### spp_gp_new assigned to the representative rank
  repr_spp <- df %>%
    inner_join(repres_df, by = c('taxon', 'spp_gp' = 'species')) %>%
    inner_join(spp_all_wide, by = c('spp_gp' = 'species')) %>%
    mutate(spp_gp_new = case_when(repres == 'genus'  ~ genus,
                                  repres == 'family' ~ family,
                                  repres == 'order'  ~ order,
                                  repres == 'class'  ~ class,
                                  TRUE ~ NA_character_)) %>%
    select(taxon, spp_gp_orig = spp_gp, spp_gp = spp_gp_new,
           everything()) %>%
    select(-phylum, -class, -order, -family, -genus) %>%
    distinct()
  
  df_higher <- df %>% 
    # filter(!spp_gp %in% repres_df$species) %>%
    left_join(rank_names_long, by = c('spp_gp' = 'name')) %>%
    rename(repres = rank)
  
  df_out <- bind_rows(df_higher, repr_spp) %>%
    select(-spp_gp_orig) %>%
    distinct()
  
  return(df_out)
}

downfill <- function(df) {
  
  if(!'repres' %in% names(df)) {
    df_repres <- assign_repr_level(df)
  } else {
    df_repres <- df
  }
  spp_all_wide <- assemble_worms(aspect = 'wide')
  
  ranks <- c('species', 'genus', 'family', 'order', 'class')
  rank_list <- vector('list', length = length(ranks)) %>%
    setNames(ranks)
  
  for(i in seq_along(ranks)) { ### i <- 4
    x <- df_repres %>%
      filter(repres == ranks[i])
    y <- x %>%
      inner_join(spp_all_wide, by = c('spp_gp' = ranks[i])) %>%
      mutate(!!ranks[i] := .data[['spp_gp']])
    rank_list[[ranks[i]]] <- y
  }
  rank_df <- bind_rows(rank_list) %>%
    distinct() %>%
    rowwise() %>%
    mutate(rank_num = which(repres == ranks)) %>%
    group_by(species) %>%
    filter(n() == 1 | rank_num == min(rank_num)) %>%
    ### this drops instances where info for a spp given at spp level, and
    ###  other info downfilled from above that may conflict w/spp level info
    ungroup() %>%
    select(-rank_num)
  
  return(rank_df)
}

##################################################.
#### function for upstream/downstream gapfill ####
##################################################.
gapfill_up_down <- function(df, col) {
  ### gf_level: 
  ### 1 = match/species, 2 = genus, 3 = family, 4 = order, 5 = class
  
  col_sd   <- paste(col, 'sd', sep = '_')
  col_nspp <- paste(col, 'nspp', sep = '_')
  
  if(!'kingdom' %in% names(df)) df$kingdom <- NA_character_
  if(!'phylum' %in% names(df)) df$phylum <- NA_character_
  if(!'db' %in% names(df)) df$db <- 'worms'
  if(!'spec_code' %in% names(df)) df$spec_code <- NA_real_
  df_tmp <- df %>%
    mutate(tmp = .data[[col]])
  gf_genus <- df_tmp %>%
    group_by(kingdom, phylum, class, order, family, genus) %>%
    summarize(across(matches('tmp'),
                     .fns = list(g_mean = ~mean(.x, na.rm = TRUE),
                                 g_sd   = ~sd(.x, na.rm = TRUE),
                                 g_cov  = ~sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE),
                                 g_nspp = ~sum(!is.na(.x))),
                     na.rm = TRUE, .names = '{.fn}_{.col}'),
              .groups = 'drop')
  gf_family <- df_tmp %>%
    group_by(kingdom, phylum, class, order, family) %>%
    summarize(across(matches('tmp'),
                     .fns = list(f_mean = ~mean(.x, na.rm = TRUE),
                                 f_sd   = ~sd(.x, na.rm = TRUE),
                                 f_cov  = ~sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE),
                                 f_nspp = ~sum(!is.na(.x))),
                     na.rm = TRUE, .names = '{.fn}_{.col}'),
              .groups = 'drop')
  gf_order <- df_tmp %>%
    group_by(kingdom, phylum, class, order) %>%
    summarize(across(matches('tmp'),
                     .fns = list(o_mean = ~mean(.x, na.rm = TRUE),
                                 o_sd   = ~sd(.x, na.rm = TRUE),
                                 o_cov  = ~sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE),
                                 o_nspp = ~sum(!is.na(.x))),
                     na.rm = TRUE, .names = '{.fn}_{.col}'),
              .groups = 'drop')
  gf_class <- df_tmp %>%
    group_by(kingdom, phylum, class) %>%
    summarize(across(matches('tmp'),
                     .fns = list(c_mean = ~mean(.x, na.rm = TRUE),
                                 c_sd   = ~sd(.x, na.rm = TRUE),
                                 c_cov  = ~sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE),
                                 c_nspp = ~sum(!is.na(.x))),
                     na.rm = TRUE, .names = '{.fn}_{.col}'),
              .groups = 'drop')
  
  df_gf_all <- df_tmp %>%
    left_join(gf_genus) %>%
    left_join(gf_family) %>%
    left_join(gf_order) %>%
    left_join(gf_class) %>%
    dplyr::select(-all_of(col)) %>%
    mutate(gf_level = case_when(!is.na(tmp) ~ 1,
                                !is.na(g_mean_tmp) ~ 2,
                                !is.na(f_mean_tmp) ~ 3,
                                !is.na(o_mean_tmp) ~ 4,
                                !is.na(c_mean_tmp) ~ 5,
                                TRUE ~ NA_real_)) %>%
    mutate(value = case_when(gf_level == 1  ~ tmp,
                             gf_level == 2  ~ g_mean_tmp,
                             gf_level == 3 ~ f_mean_tmp,
                             gf_level == 4  ~ o_mean_tmp,
                             gf_level == 5  ~ c_mean_tmp,
                             TRUE ~ NA_real_)) %>%
    mutate(sd = case_when(gf_level == 1  ~ NA_real_,
                          gf_level == 2  ~ g_sd_tmp,
                          gf_level == 3 ~ f_sd_tmp,
                          gf_level == 4  ~ o_sd_tmp,
                          gf_level == 5  ~ c_sd_tmp,
                          TRUE ~ NA_real_)) %>%
    mutate(nspp = case_when(gf_level == 1  ~ NA_integer_,
                            gf_level == 2  ~ g_nspp_tmp,
                            gf_level == 3 ~ f_nspp_tmp,
                            gf_level == 4  ~ o_nspp_tmp,
                            gf_level == 5  ~ c_nspp_tmp,
                            TRUE ~ NA_integer_)) %>%
    ungroup() %>%
    mutate(trait = col) %>%
    dplyr::select(db, spec_code, 
                  kingdom, phylum, class, order, family, genus, species, 
                  trait, value, sd, nspp, gf_level) %>%
    distinct()
  
  return(df_gf_all)
}

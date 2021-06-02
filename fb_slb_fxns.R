###################################
### Helper functions in general ###
###################################
here_anx <- function(f) { 
  ### create file path to git-annex dir for project
  f_anx <- sprintf('~/git-annex/spp_vuln_mapping/%s', f)
}

#####################################################
### Helper functions for FishBase and SeaLifeBase ###
#####################################################
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
      select(spec_code, species, db, keep_fxn(keep_cols))
    slb <- slb %>%
      select(spec_code, species, db, keep_fxn(keep_cols))
  }
  if(!is.null(drop_cols)) {
    fb <- fb %>%
      select(-drop_fxn(drop_cols))
    slb <- slb %>%
      select(-drop_fxn(drop_cols))
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

#################################################
### function for downstream direct match fill ###
#################################################

downfill <- function(df) {
  x <- data.table::fread(here_anx('vuln_gapfilled_all_tx.csv')) %>%
    select(-vuln_gf_id) %>%
    filter_dupe_spp(level = 'class')
  
  ranks <- c('species', 'genus', 'family', 'order', 'class')
  
  rank_list <- vector('list', length = length(ranks)) %>%
    setNames(ranks)
  for(rank in ranks) { ### rank <- ranks[2]
    xx <- x %>%
      filter(match_rank == rank)
    y <- df %>%
      right_join(xx, by = c('spp_gp' = rank, 'taxon')) %>%
      mutate(!!rank := .data[['spp_gp']])
    rank_list[[rank]] <- y
  }
  rank_df <- bind_rows(rank_list) %>%
    select(taxon, spp_gp, # category, 
           class, order, family, genus, species, 
           gapfill, match_rank, everything())
  
  return(rank_df)
}

################################################
### function for upstream/downstream gapfill ###
################################################
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
    select(-all_of(col)) %>%
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
    select(db, spec_code, 
           kingdom, phylum, class, order, family, genus, species, 
           trait, value, sd, nspp, gf_level) %>%
    distinct()
  
  return(df_gf_all)
}


###################################################################
### assemble overall species taxonomy from extracted WoRMS data ###
###################################################################
get_worms <- function(animalia_only = TRUE) {
  worms_f <- here_anx('worms_taxa_animalia.csv')
  if(!file.exists(worms_f)) {
    
    worms_p_from_k <- read_csv(here('~/git-annex/spp_vuln_mapping/taxa_from_worms/phylum_from_kingdom_worms.csv')) %>%
      select(phylum = name, kingdom = parent)
    worms_c_from_p <- read_csv(here('~/git-annex/spp_vuln_mapping/taxa_from_worms/class_from_phylum_worms.csv')) %>%
      select(class = name, phylum = parent)
    worms_o_from_c <- read_csv(here('~/git-annex/spp_vuln_mapping/taxa_from_worms/order_from_class_worms.csv')) %>%
      select(order = name, class = parent)
    worms_f_from_o <- read_csv(here('~/git-annex/spp_vuln_mapping/taxa_from_worms/family_from_order_worms.csv')) %>%
      select(family = name, order = parent)
    worms_g_from_f <- read_csv(here('~/git-annex/spp_vuln_mapping/taxa_from_worms/genus_from_family_worms.csv')) %>%
      select(genus = name, family = parent)
    worms_s_from_g <- read_csv(here('~/git-annex/spp_vuln_mapping/taxa_from_worms/species_from_genus_worms.csv')) %>%
      select(species = name, genus = parent)
    
    spp_from_worms_raw <- worms_p_from_k %>%
      full_join(worms_c_from_p) %>%
      full_join(worms_o_from_c) %>%
      full_join(worms_f_from_o) %>%
      full_join(worms_g_from_f) %>%
      full_join(worms_s_from_g) %>%
      filter(!is.na(species)) %>%
      select(kingdom, phylum, class, order, family, genus, species) 
    
    ### manually fix incorrect gapfilled taxa
    ### filter: left side is the troublesome lower-level taxon, 
    ### right side is the accurate upper-level taxon
    spp_from_worms <- spp_from_worms_raw %>%
      filter_dupe_spp()
    # x <- check_dupe_spp(spp_from_worms)
    
    write_csv(spp_from_worms, worms_f)
    return(spp_from_worms)
  } else {
    return(read_csv(worms_f))
  }
}
  

check_dupe_spp <- function(df) {
  x <- show_dupes(df, 'species') %>%
    left_join(read_csv(here('_data/spp_traits_valid.csv')) %>% 
                downfill() %>%
                select(taxon, species) %>%
                distinct()) %>%
    group_by(species) %>%
    summarize(n_phyla = n_distinct(phylum),
              phyla   = paste(unique(phylum), collapse = ';'),
              n_class = n_distinct(class),
              classes = paste(unique(class), collapse = ';'),
              orders = paste(unique(order), collapse = ';'),
              families = paste(unique(family), collapse = ';'),
              genus = paste(unique(genus), collapse = ';'))
  return(x)
}

filter_dupe_spp <- function(df, level = c('all', 'class')[1]) {
  df_out <- df %>%
    filter(!(genus %in% c('aturia', 'pleurifera', 'distefanella', 'philippiella') & class != 'mollusca')) %>%
    filter(!(family == 'margaritidae' & class  != 'gastropoda')) %>%
    filter(!(genus == 'spongicola' & class != 'malacostraca')) %>%
    filter(!(genus == 'trachyaster' & class != 'echinoidea')) %>%
    filter(!(genus == 'nematoporella' & class != 'stenolaemata')) %>%
    filter(!(genus == 'tremaster' & class != 'asteroidea')) %>%
    filter(!(genus == 'versuriga' & family != 'versurigidae')) %>%
    filter(!(genus == 'diplocoenia' & family != 'faviidae')) %>%
    filter(!(genus == 'stictostega' & family != 'hippothoidae')) %>%
    filter(!(genus == 'cladochonus' & family != 'auloporidae')) %>%
    distinct()
  
  if(level == 'all') {
    df_out <- df_out %>%
      filter(kingdom == 'animalia') %>%
      filter(!(genus  == 'eugyra' & phylum != 'chordata')) %>%
      filter(!(genus %in% c('bathya', 'thoe') & phylum != 'arthropoda')) %>%
      filter(!(genus %in% c('bergia', 'geminella') & phylum != 'cnidaria')) %>%
      filter(!(genus %in% c('chaperia', 'favosipora') & phylum != 'bryozoa')) %>%
      filter(!(genus == 'atractotrema' & phylum != 'platyhelminthes')) %>%
      filter(!(genus == 'ctenella' & phylum != 'ctenophora'))
  }
  return(df_out)
}


### Get vulnerability traits by species
get_vuln_traits <- function() {
  read_csv(here('_data/traits_vulnerability/spp_traits_valid.csv'))
}

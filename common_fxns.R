###################################
### Helper functions in general ###
###################################
here_anx <- function(f = '', ...) { 
  ### create file path to git-annex dir for project
  f <- paste(f, ..., sep = '/')
  f <- stringr::str_replace_all(f, '\\/+', '/')
  f_anx <- sprintf('/home/shares/ohi/spp_vuln/spp_vuln_mapping/%s', f)
  return(f_anx)
}

#####################################
### Helper functions for AquaMaps ###
#####################################

am_dir <- '/home/shares/ohi/spp_vuln/aquamaps_2021'
get_am_spp_info <- function() {
  ### resolves AquaMaps species names according to WoRMS classification,
  ### drops subspp names
  am_spp_resolve <- data.table::fread(here('_data/worms_taxa/aquamaps_aphia_records.csv')) %>%
    filter(!is.na(aphia_id)) %>%
    group_by(am_sciname) %>%
    mutate(match = am_sciname == valid_name) %>%
    arrange(desc(match)) %>% ### get TRUE on top, if any
    summarize(valid_name = first(valid_name), aphia_id = first(aphia_id))
  
  spp_info <- data.table::fread(file.path(am_dir, 'ver10_2019_speciesoccursum_iucn.csv')) %>%
    janitor::clean_names() %>%
    rename(am_sid = species_id, iucn_sid = iucn_id, comname = f_bname) %>%
    mutate(am_sciname = tolower(paste(genus, species))) %>%
    full_join(am_spp_resolve, by = 'am_sciname') %>%
    mutate(sciname = ifelse(is.na(valid_name), am_sciname, valid_name)) %>%
    mutate(sciname = case_when(str_detect(sciname, 'incertae sedis') ~ sciname,
                               str_detect(sciname, ' var\\.') ~ str_remove(sciname, ' var\\. .+'),
                               str_detect(sciname, '\\(|\\)') ~ str_extract(sciname, '[^ ]+ \\(.+?\\) [^ ]+'),
                               TRUE ~ str_extract(sciname, '[^ ]+ [^ ]+'))) %>%
    mutate(across(where(is.character), ~tolower(.x))) %>%
    clean_scinames('sciname')
  
  return(spp_info)
}

get_am_spp_envelopes <- function() {
  env_file <- file.path(am_dir, 'species_envelope_summary.csv')
  
  ### create a tidy summary of environmental envelopes
  # if(!file.exists) { ### create from scratch
  #   env_info <- data.table::fread(file.path(am_dir, 'ver10_2019_hspen.csv')) %>%
  #     janitor::clean_names() %>%
  #     rename(am_sid = species_id)
  #   
  #   params <- c('depth', 'temp', 'salinity', 'prim_prod', 'ice_con', 'oxy', 'land_dist')
  #   prm_or <- paste0(params, collapse = '|')
  #   dist_vals <- c('min', 'pref_min', 'mean', 'pref_max', 'max')
  #   env_sum_df <- env_info %>%
  #     rename(depth_mean = mean_depth) %>%
  #     select(am_sid, contains(params)) %>%
  #     gather(param_full, value, -am_sid) %>%
  #     mutate(param = str_extract(param_full, prm_or),
  #            dist = str_replace(param_full, paste0('(', prm_or, ')_'), '')) %>%
  #     select(-param_full) %>%
  #     spread(dist, value) %>%
  #     filter(yn == 1) %>%
  #     select(-yn) %>%
  #     gather(dist, value, -am_sid, -param) %>%
  #     mutate(dist = factor(dist, levels = dist_vals),
  #            value = ifelse(param == 'depth', log10(value + 1), value),
  #            param = ifelse(param == 'depth', 'log10_depth', param))
  #   
  #   write_csv(env_sum_df, env_file)
  # }
  return(data.table::fread(env_file))
}

get_hcaf_info <- function() {
  hcaf_info <- data.table::fread(file.path(am_dir, 'hcaf_v7.csv')) %>%
    janitor::clean_names() %>%
    mutate(across(where(is.numeric), ~ ifelse(.x == -9999, NA, .x)))
  return(hcaf_info)
}

get_hcaf_rast <- function() {
  rast_base <- raster::raster(ext = raster::extent(c(-180, 180, -90, 90)), 
                              resolution = 0.5, crs = '+init=epsg:4326') %>%
    raster::setValues(1:raster::ncell(.))
  return(rast_base)
}

map_to_hcaf <- function(df, by = 'loiczid', which, xfm = NULL) {
  if(!by %in% names(df)) stop('Dataframe needs a valid column for "by" (e.g., loiczid)!')
  if(!which %in% names(df)) stop('Dataframe needs a valid column for "which" (e.g., n_spp)!')
  if(any(is.na(df[[by]]))) stop('Dataframe contains NA values for "by"!')
  
  r_base <- get_hcaf_rast()
  
  out_rast <- raster::subs(x = r_base, y = df, 
                           by = by, which = which)
  
  if(!is.null(xfm)) {
    if(class(xfm) != 'function') stop('xfm argument must be a function!')
    out_rast <- xfm(out_rast)
  }
  return(out_rast)
}

get_am_spp_cells <- function(occurcells_cut = 10, prob_cut = 0) {
  
  spp_cell_file <- file.path(am_dir, 'hcaf_species_native_clean.csv')
  
  ### create a cleaner version  of spp native file for speed and size!
  # if(!file.exists(spp_cell_file)) {
  #   csq_loiczid <- get_hcaf_info() %>%
  #     select(loiczid, csquare_code) %>%
  #     distinct()
  #   spp_cells <- data.table::fread(file.path(am_dir, 'hcaf_species_native.csv')) %>%
  #     janitor::clean_names() %>%
  #     oharac::dt_join(csq_loiczid, by = 'csquare_code', type = 'left') %>%
  #     select(am_sid = species_id, loiczid, prob = probability) %>%
  #     distinct() %>%
  #     mutate(am_sid = tolower(am_sid))
  # 
  #   write_csv(spp_cells, spp_cell_file)
  # }
  
  spp_cells <- data.table::fread(spp_cell_file)
  
  if(occurcells_cut > 0) {
    spp_valid <- get_am_spp_info() %>%
      filter(occur_cells >= occurcells_cut) %>%
      dplyr::select(am_sid, sciname)
    spp_cells <- spp_cells %>%
      filter(am_sid %in% spp_valid$am_sid)
  }
  if(prob_cut > 0) {
    spp_cells <- spp_cells %>%
      filter(prob >=  prob_cut)
  }
  return(spp_cells)
}

#####################################################
###    Helper functions for vulnerability data    ###
#####################################################
get_spp_vuln <- function(gapfill = c('family', 'all')[2]) {
  if(gapfill == 'family') {
    ### read from pre-filtered gapfill files in Github
    vuln_tx <- data.table::fread(here('_data/vuln_data/vuln_gapfilled_tx.csv'))
    vuln_score <- data.table::fread(here('_data/vuln_data/vuln_gapfilled_score.csv')) %>%
      gather(stressor, score, -vuln_gf_id)
    vuln_sd <- data.table::fread(here('_data/vuln_data/vuln_gapfilled_sd.csv')) %>%
      gather(stressor, sd_score, -vuln_gf_id)
  } else {
    ### read from unfiltered gapfilled files on Mazu server
    mazu_dir <- '/home/shares/ohi/spp_vuln/spp_vuln_framework/post_gapfill'
    vuln_tx <- data.table::fread(file.path(mazu_dir, 'vuln_gapfilled_all_tx.csv'))
    vuln_score <- data.table::fread(file.path(mazu_dir, 'vuln_gapfilled_all_score.csv')) %>%
      gather(stressor, score, -vuln_gf_id)
    vuln_sd <- data.table::fread(file.path(mazu_dir, 'vuln_gapfilled_all_sd.csv')) %>%
      gather(stressor, sd_score, -vuln_gf_id)
  }
  vuln_df <- vuln_tx %>%
    oharac::dt_join(vuln_score, by = c('vuln_gf_id'), type = 'full') %>%
    oharac::dt_join(vuln_sd,    by = c('vuln_gf_id', 'stressor'), type = 'full') %>%
    dplyr::select(-vuln_gf_id) %>%
    clean_scinames('species') %>%
    clean_scinames('genus') %>%
    filter(!is.na(score))
  return(vuln_df)
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

#################################################
### function for downstream direct match fill ###
#################################################

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


###################################################################
### assemble overall species taxonomy from extracted WoRMS data ###
###################################################################
assemble_worms <- function(aspect = 'wide', seabirds_only = TRUE, am_patch = TRUE) {
  ### Note: this drops all kingdoms but Animalia 
  
  p_from_k <- data.table::fread(here('_data/worms_taxa', 
                            'expand1_phylum_from_kingdom_worms.csv')) %>%
    filter(!is.na(id)) %>%
    select(-id) %>%
    distinct() %>%
    mutate
  c_from_p <- data.table::fread(here('_data/worms_taxa', 
                            'expand2_class_from_phylum_worms.csv')) %>%
    filter(!is.na(id)) %>%
    select(-id) %>%
    distinct()
  o_from_c <- data.table::fread(here('_data/worms_taxa', 
                            'expand3_order_from_class_worms.csv')) %>%
    filter(!is.na(id)) %>%
    select(-id) %>%
    distinct()
  f_from_o <- data.table::fread(here('_data/worms_taxa', 
                            'expand4_family_from_order_worms.csv')) %>%
    filter(!is.na(id)) %>%
    select(-id) %>%
    distinct()
  g_from_f <- data.table::fread(here('_data/worms_taxa', 
                            'expand5_genus_from_family_worms.csv')) %>%
    filter(!is.na(id)) %>%
    select(-id) %>%
    distinct()
  s_from_g <- data.table::fread(here('_data/worms_taxa', 
                            'expand6_species_from_genus_worms.csv')) %>%
    filter(!is.na(id)) %>%
    select(-id) %>%
    distinct()
  
  if(am_patch) {
    am_patch_wide <- data.table::fread(here('_data/worms_taxa',
                                   'expand7_aquamaps_patch.csv')) %>%
      distinct() %>% 
      mutate(source = 'am')
  } else {
    am_patch_wide <- data.frame(source = 'am') ### blank dataframe for bind_rows
  }
  
  rank_lvls <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
  
  ### create wide for complete classification for each species
  spp_wide <- s_from_g %>%
    select(genus = parent, species = name) %>%
    left_join(g_from_f %>% select(family = parent, genus = name), 
              by = c('genus')) %>%
    left_join(f_from_o %>% select(order = parent, family = name), 
              by = c('family')) %>%
    left_join(o_from_c %>% select(class = parent, order = name), 
              by = c('order')) %>%
    left_join(c_from_p %>% select(phylum = parent, class = name), 
              by = c('class')) %>%
    left_join(p_from_k %>% select(kingdom = parent, phylum = name),
              by = c('phylum')) %>%
    select(phylum, class, order, family, genus, species) %>%
    mutate(source = 'worms') %>%
    bind_rows(am_patch_wide) %>%
    filter(phylum %in% p_from_k$name) %>%
    ### since p_from_k only includes Animalia, this phylum selection drops
    ### AquaMaps non-animalia phyla, e.g., chlorophytes, cyanobacteria, plants
    clean_scinames('species') %>%
    clean_scinames('genus') %>%
    clean_scinames('family') %>%
    clean_scinames('order') %>%
    clean_scinames('class') %>%
    clean_scinames('phylum') %>%
    distinct()
  
  spp_wide_disambiguated <- disambiguate_species(spp_wide)
  
  spp_wide_am_resolved <- resolve_am_disputes(spp_wide_disambiguated)
  
  ### disambiguate "not assigned" - these appear at order and class levels; 
  ### this would only matter when doing upstream/downstream imputation at 
  ### order/class levels...
  spp_df <- spp_wide_am_resolved %>%
    mutate(class = ifelse(class == 'not assigned', paste(phylum, 'not assigned'), class),
           order = ifelse(order == 'not assigned', paste(class, 'not assigned'), order))
  
  
  if(seabirds_only == TRUE) {
    # seabird_list <- readxl::read_excel(here('_raw/species_numbers.xlsx'),
    #                                    sheet = 'seabirds', skip = 1) %>%
    #   janitor::clean_names() %>%
    #   select(spp = scientific_name_fixed) %>%
    #   filter(!is.na(spp)) %>%
    #   mutate(spp = tolower(spp) %>% str_trim) %>%
    #   select(spp) %>%
    #   distinct()
    seabird_list <- read_csv(here('_data/seabird_list.csv'), show_col_types = FALSE)
    spp_df <- spp_df %>%
      filter(!(tolower(class) == 'aves' & !tolower(species) %in% seabird_list$spp))
  }
  
  if(aspect == 'long') {
    
    ### make that long, but keeping structure for each species
    spp_df_long <- spp_df %>%
      mutate(spp = species) %>%
      gather(rank, name, phylum:species) %>%
      mutate(rank = factor(rank, levels = rank_lvls))
    return(spp_df_long)
    
  } else {
    
    return(spp_df)
    
  }
}

resolve_am_disputes <- function(spp_wide) {
  ## coerce AM classifications to match WoRMS
  dupes_g <- show_dupes(spp_wide %>%
                          select(-species) %>%
                          distinct(),
                        'genus') %>%
    filter(!(genus == 'gammarus' & family == 'paratanaidae')) %>%
    group_by(genus) %>%
    filter(n_distinct(class) > 1 | n_distinct(order) > 1 | n_distinct(family) > 1) %>%
    mutate(am_count = sum(source == 'am'),
           non_am_count = sum(source == 'worms'))
  ### identify duplicated genera that are included (uniquely) in WoRMS
  worms_fill <- dupes_g %>%
    filter(non_am_count > 0) %>%
    filter(source == 'worms')
  
  dupes_g_force_worms <- spp_wide %>%
    filter(genus %in% worms_fill$genus) %>%
    rowwise() %>%
    mutate(family = worms_fill$family[genus == worms_fill$genus],
           order  = worms_fill$order[genus == worms_fill$genus],
           class  = worms_fill$class[genus == worms_fill$genus]) %>%
    ungroup()
  
  nonworms_fill <- dupes_g %>%
    filter(!genus %in% dupes_g_force_worms$genus) %>%
    mutate(genus  = case_when(genus == 'lamellaria'    ~ 'lamellidea',
                              genus == 'leptocephalus' ~ 'conger',
                              TRUE ~ genus)) %>%
    mutate(family = case_when(genus == 'conger'        ~ 'congridae',
                              genus == 'cerithiella'   ~ 'newtoniellidae',
                              genus == 'elysia'        ~ 'plakobranchidae',
                              genus == 'euciroa'       ~ 'euciroidae',
                              genus == 'eulimastoma'   ~ 'pyramidellidae',
                              genus == 'lamellidea'    ~ 'achatinellidae',
                              genus == 'mathilda'      ~ 'mathildidae',
                              genus == 'polybranchia'  ~ 'hermaeidae',
                              genus == 'tjaernoeia'    ~ 'tjaernoeiidae',
                              genus == 'lamellitrochus' ~ 'solariellidae',
                              TRUE ~ family)) %>%
    mutate(order = case_when(genus == 'lamellitrochus' ~ 'trochida',
                             genus == 'elysia' ~ 'sacoglossa',
                             genus == 'eulimastoma' ~ 'pylopulmonata',
                             genus == 'conger' ~ 'anguilliformes',
                             genus == 'polybranchia' ~ 'sacoglossa',
                             genus == 'tjaernoeia' ~ '[unassigned] euthyneura',
                             TRUE ~ order)) %>%
    mutate(class = case_when(order == 'sacoglossa' ~ 'gastropoda',
                             TRUE ~ class)) %>%
    mutate(phylum = case_when(class == 'gastropoda' ~ 'mollusca',
                             TRUE ~ phylum)) %>%
    distinct()
  
  dupes_g_force_nonworms <- spp_wide %>%
    filter(genus %in% dupes_g$genus & !genus %in% worms_fill$genus) %>%
    rowwise() %>%
    mutate(genus  = case_when(genus == 'lamellaria'    ~ 'lamellidea',
                              genus == 'leptocephalus' ~ 'conger',
                              TRUE ~ genus)) %>%
    mutate(family = nonworms_fill$family[genus == nonworms_fill$genus],
           order  = nonworms_fill$order[genus  == nonworms_fill$genus],
           class  = nonworms_fill$class[genus  == nonworms_fill$genus],
           phylum = nonworms_fill$phylum[genus == nonworms_fill$genus]) %>%
    ungroup() %>%
    distinct()
  
  spp_wide_am_resolved <- spp_wide %>%
    filter(!genus %in% dupes_g$genus) %>%
    bind_rows(dupes_g_force_worms, dupes_g_force_nonworms) %>%
    select(-source) %>%
    distinct()
  
  return(spp_wide_am_resolved)
}

disambiguate_species <- function(spp_wide) {
  
  dupes <- spp_wide %>%
    oharac::show_dupes('species')
  dupes_drop_source <- dupes %>%
    select(-source) %>%
    distinct() %>%
    show_dupes('species')
  ### currently, none due to appearing in both AM and WoRMS
  
  spp_wide_nodupes <- spp_wide %>%
    filter(!species %in% dupes_drop_source$species)
  
  dupes_fixed <- dupes %>%
    filter(species != 'no match') %>%
    mutate(genus  = ifelse(species == 'praephiline finmarchica', 'praephiline', genus),
           family = ifelse(species == 'praephiline finmarchica', 'laonidae', family),
           genus  = ifelse(species == 'polititapes rhomboides', 'polititapes', genus)) %>%
    mutate(keep = case_when(family == 'margaritidae' & order == 'trochida' & genus != 'pinctada' ~ TRUE,
                            genus  == 'pinctada'     & order == 'ostreida'       ~ TRUE,
                            genus  == 'atractotrema' & class == 'gastropoda'     ~ TRUE,
                            genus  == 'chaperia'     & phylum == 'bryozoa'       ~ TRUE,
                            family == 'molgulidae'   & genus == 'eugyra'         ~ TRUE,
                            genus  == 'aturia'       & order == 'nautilida'      ~ TRUE,
                            genus  == 'spongicola'   & family == 'spongicolidae' ~ TRUE,
                            genus  == 'stictostega'  & family == 'hippothoidae'  ~ TRUE,
                            genus  == 'favosipora'   & family == 'densiporidae'  ~ TRUE,
                            genus  == 'cladochonus'  & family == 'pyrgiidae'     ~ TRUE,
                            genus  == 'bathya'       & order == 'amphipoda'      ~ TRUE,
                            genus  == 'ctenella'     & family == 'ctenellidae'   ~ TRUE,
                            genus  == 'pleurifera'   & family == 'columbellidae' ~ TRUE,
                            genus  == 'thoe'         & family == 'mithracidae'   ~ TRUE,
                            genus  == 'diplocoenia'  & family == 'acroporidae'   ~ TRUE,
                            genus  == 'versuriga'    & family == 'versurigidae'  ~ TRUE,
                            genus  == 'tremaster'    & family == 'asterinidae'   ~ TRUE,
                            genus  == 'distefanella' & family == 'radiolitidae'  ~ TRUE,
                            genus  == 'bracthelia'   & family == 'agatheliidae'  ~ TRUE,
                            genus  == 'dinetia'      & family == 'draconematidae'   ~ TRUE,
                            genus  == 'bergia'       & family == 'drepanophoridae'  ~ TRUE,
                            genus  == 'geminella'    & family == 'catenicellidae'   ~ TRUE,
                            genus  == 'nematoporella' & family == 'arthrostylidae'  ~ TRUE,
                            genus  == 'philippiella' & family == 'steinmanellidae'  ~ TRUE,
                            genus  == 'trachyaster'  & family == 'palaeostomatidae' ~ TRUE,
                            TRUE ~ FALSE)) %>%
    group_by(species) %>%
    mutate(gen_match = str_detect(species, paste0('^', genus, ' ')) & source == 'am',
           keep = ifelse(sum(!keep) > 1 & gen_match, TRUE, keep)) %>%
    filter(keep) %>%
    select(-keep, -gen_match) %>%
    distinct()

  spp_wide_clean <- bind_rows(spp_wide_nodupes, dupes_fixed)
  
  return(spp_wide_clean)
}


check_dupe_spp <- function(df) {
  x <- show_dupes(df, 'species') %>%
    left_join(data.table::fread(here('_data/spp_traits_valid.csv')) %>% 
                downfill() %>%
                dplyr::select(taxon, species) %>%
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
  x <- data.table::fread(here('_data/traits_vulnerability/spp_traits_valid.csv'))
}

####################################################=
####  Resolve names using WoRMS API fuzzy match  ####
####################################################=

fuzzy_match <- function(spp_vec, marine_only = TRUE) {
  
  spp_fix <- str_replace_all(spp_vec, ' +', '%20')
  spp_arg <- paste0('scientificnames[]=', spp_fix, collapse = '&')
  
  mar_flag <- tolower(as.character(marine_only))
  matchname_stem <- 'https://www.marinespecies.org/rest/AphiaRecordsByMatchNames?%s&marine_only=%s'
  matchname_url  <- sprintf(matchname_stem, spp_arg, mar_flag)
  
  Sys.sleep(0.25) ### slight pause for API etiquette
  match_records <- try(jsonlite::fromJSON(matchname_url)) 
  if(class(match_records) == 'try-error') {
    match_df <- data.frame(orig_sciname = spp_vec,
                           valid_name = 'no match',
                           aphia_id = -9999) 
    return(match_df)
  } else {
    match_df <- match_records %>%
      setNames(spp_vec) %>%
      bind_rows(.id = 'orig_sciname')
    
    still_unmatched <- data.frame(orig_sciname = spp_vec) %>%
      filter(!orig_sciname %in% match_df$orig_sciname) %>%
      mutate(valid_name = 'no match', aphia_id = -9999)
    out_df <- match_df %>%
      bind_rows(still_unmatched) %>%
      select(orig_sciname, valid_name, aphia_id = valid_AphiaID) %>%
      mutate(valid_name = tolower(valid_name))
    return(out_df)
  }
}

collect_records <- function(fb_df, field, file_tag, force = FALSE) {
  
  ### get names with values that don't match WoRMS names
  spp_from_worms <- assemble_worms()
  
  fb_df <- fb_df %>%
    rename(tmp := !!field)
  
  no_match <- anti_join(fb_df, spp_from_worms, by = 'species')
  no_match_w_vals <- no_match %>%
    filter(!is.na(tmp)) %>%
    mutate(genus = str_extract(species, '^[a-z]+(?= ?)'),
           genus_match = genus %in% spp_from_worms$genus)
  
  # table(no_match_w_vals %>% select(genus_match, db))
  
  names_to_resolve <- no_match_w_vals %>%
    filter(genus_match) %>%
    .$species %>% unique() %>% sort()
  
  ### Define file and check whether it exists
  aphia_records_csv <- sprintf(here('int/%s.csv'), file_tag)
  
  
  if(!file.exists(aphia_records_csv) | force) {
    
    chunk_size <- 25
    n_chunks <- ceiling(length(names_to_resolve) / chunk_size)
    record_chunk_stem <- 'tmp/%s_chunk_%s_%s.csv'
    
    if(force) {
      unlink(list.files(here('tmp'), pattern = sprintf('%s_chunk', file_tag), full.names = TRUE))
    }
    
    for(i in 1:n_chunks) { ### i <- 1
      message('Processing chunk ', i, ' of ', n_chunks)
      i_start <- (i-1) * chunk_size + 1
      i_end   <- min(i * chunk_size, length(names_to_resolve))
      chunk_csv <- here(sprintf(record_chunk_stem, file_tag, i_start, i_end))
      if(file.exists(chunk_csv)) {
        message('Chunk exists: ', basename(chunk_csv), '... skipping!')
        next()
      }
      spp_vec <- names_to_resolve[i_start:i_end]
      chunk_df <- fuzzy_match(spp_vec, marine_only = FALSE)
      write_csv(chunk_df, chunk_csv)
    }
    
    chunk_fs <- list.files(here('tmp'), pattern = sprintf('%s_chunk', file_tag), full.name = TRUE)
    record_df <- parallel::mclapply(chunk_fs, data.table::fread) %>%
      bind_rows() %>%
      distinct()
    write_csv(record_df, aphia_records_csv)
  }
  record_df <- data.table::fread(aphia_records_csv) %>%
    clean_scinames('valid_name')
}

clean_scinames <- function(df, field) {
  ### eliminate confounding clutter in taxonomic names
  df_clean <- df %>%
    rename(tmp := !!field) %>%
    mutate(tmp = str_remove_all(tmp, '\\(.+?\\)'), ### parentheticals
           tmp = str_remove_all(tmp, '\\[.+?\\]'), ### brackets
           tmp = str_remove_all(tmp, '[^a-z ]')) %>% ### punctuation
    mutate(tmp = str_squish(tmp)) %>%
    rename(!!field := tmp)
}

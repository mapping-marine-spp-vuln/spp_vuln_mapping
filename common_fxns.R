### set options
options(readr.show_types = FALSE)

#####################################.
#### Helper functions in general ####
#####################################.
here_anx <- function(f = '', ...) { 
  ### create file path to git-annex dir for project
  f <- paste(f, ..., sep = '/')
  f <- stringr::str_replace_all(f, '\\/+', '/')
  f_anx <- sprintf('/home/shares/ohi/spp_vuln/spp_vuln_mapping/%s', f)
  return(f_anx)
}

smash2or <- function(x, sep = ' ') {
  ### collapse to OR string
  paste(x, sep = sep, collapse = '|')
}

assemble_spp_info_df <- function(fe_only = TRUE, vuln_only = TRUE) {
  spp_am <- get_am_spp_info()  %>%
    filter(occur_cells >= 10) %>%
    select(species = sciname) %>%
    mutate(am_mapped = TRUE) %>%
    distinct()
  
  spp_iucn <- read_csv(here('_data/iucn_spp/iucn_to_worms_match.csv'), 
                       show_col_types = FALSE) %>%
    rename(species = worms_name, iucn_mapped = mapped)
  
  spp_vuln <- get_spp_vuln() %>%
    rename(v_score = score)
  
  spp_worms <- assemble_worms() %>%
    select(species) %>%
    distinct()
  
  spp_fe <- get_fe_traits()
  
  spp_map_fstem <- here_anx('spp_maps_mol', '%s_spp_mol_%s.csv')
  
  spp_df <- spp_worms %>%
    oharac::dt_join(spp_vuln, by = 'species', type = 'left') %>%
    oharac::dt_join(spp_am,   by = 'species', type = 'left') %>%
    oharac::dt_join(spp_iucn, by = 'species', type = 'left') %>%
    oharac::dt_join(spp_fe, by = 'species', 
                    type = ifelse(fe_only, 'inner', 'left')) %>%
    filter(!is.na(v_score) | !vuln_only) %>%
    ### mapped in AquaMaps and/or IUCN
    filter(am_mapped | iucn_mapped) %>%
    mutate(src = ifelse(iucn_mapped & !is.na(iucn_mapped), 'iucn', 'am'),
           id  = ifelse(src == 'iucn', iucn_sid, str_replace_all(species, ' ', '_')),
           map_f = sprintf(spp_map_fstem, src, id)) %>%
    distinct()
  
  return(spp_df)
}

### Helper functions for gathering species rangemaps
check_tryerror <- function(l) {
  x <- sapply(l, class) %>% 
    unlist() %>% as.vector()
  return(any(stringr::str_detect(tolower(x), 'error')))
}

get_one_map <- function(f) {
  if(file.exists(f)) {
    df <- data.table::fread(f)
    
    if('presence' %in% names(df)) {   ### IUCN spp map with presence col
      df <- df %>% filter(presence != 5) %>% select(-presence)
    }
    if('prob' %in% names(df)) {
      df <- df %>% filter(prob >= 0.5) %>% select(-prob)
    }
    return(df)
  } else {
    warning('No map found for ', basename(f))
    return(NULL)
  }
}

collect_spp_rangemaps <- function(spp_vec, file_vec, idcol = 'species', parallel = TRUE) {
  ### give a vector of species names (or IDs) and filenames; 
  ### default: read in using parallel::mclapply
  message('Collecting ', n_distinct(file_vec), ' maps for ', 
          n_distinct(spp_vec), ' species...')
  if(parallel == TRUE) {
    out_maps_list <- parallel::mclapply(file_vec, mc.cores = 40, FUN = get_one_map)
  } else {
    out_maps_list <- lapply(file_vec, FUN = get_one_map)
  }
  if(check_tryerror(out_maps_list)) {
    stop('Try-error found when assembling species rangemaps!')
  }
  message('... Binding maps...')
  out_maps_df <- out_maps_list %>%
    setNames(spp_vec) %>%
    purrr::compact() %>%
    data.table::rbindlist(idcol = idcol) %>%
    distinct()
  
  return(out_maps_df)
}


#######################################.
#### Helper functions for AquaMaps ####
#######################################.

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
  env_df <- data.table::fread(env_file) %>%
    mutate(am_sid = tolower(am_sid))
  return(env_df)
}

get_hcaf_info <- function() {
  hcaf_info <- data.table::fread(file.path(am_dir, 'hcaf_v7.csv')) %>%
    janitor::clean_names() %>%
    mutate(across(where(is.numeric), ~ ifelse(.x == -9999, NA, .x)))
  return(hcaf_info)
}

get_hcaf_rast <- function() {
  rast_base <- terra::rast(ext(c(-180, 180, -90, 90)), 
                           resolution = 0.5, crs = '+init=epsg:4326') %>%
    terra::setValues(1:terra::ncell(.))
  return(rast_base)
}

map_to_hcaf <- function(df, by = 'loiczid', which, xfm = NULL, ocean_mask = FALSE) {
  if(!by %in% names(df)) stop('Dataframe needs a valid column for "by" (e.g., loiczid)!')
  if(!which %in% names(df)) stop('Dataframe needs a valid column for "which" (e.g., n_spp)!')
  if(any(is.na(df[[by]]))) stop('Dataframe contains NA values for "by"!')
  if(length(df[[by]]) != length(unique(df[[by]]))) stop('Dataframe contains duplicate observations of ', by, '!')
  
  out_rast <- get_hcaf_rast()
  
  df1 <- data.frame(cell_id = 1:ncell(out_rast)) %>%
    rename(!!by := cell_id) %>%
    dt_join(df, by = by, type = 'left')
  values(out_rast) <- df1[[which]]

  if(!is.null(xfm)) {
    if(class(xfm) != 'function') stop('xfm argument must be a function!')
    out_rast <- xfm(out_rast)
  }

  if(ocean_mask) {
    out_rast <- out_rast %>%
      terra::mask(terra::rast(here('_spatial/ocean_area_wgs84_0.5deg.tif')))
  }
  return(out_rast)
}

get_am_spp_cells <- function(occurcells_cut = 0, prob_cut = 0) {
  
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

###############################################.
####   Helper functions for IUCN species   ####
###############################################.

get_mol_rast <- function() {
  rast_base <- terra::rast(here('_spatial/ocean_area_mol.tif')) %>%
    terra::setValues(1:terra::ncell(.))
  return(rast_base)
}

map_to_mol <- function(df, by = 'cell_id', which, xfm = NULL, ocean_mask = TRUE) {
  if(!by %in% names(df)) stop('Dataframe needs a valid column for "by" (e.g., cell_id)!')
  if(!which %in% names(df)) stop('Dataframe needs a valid column for "which" (e.g., n_spp)!')
  if(any(is.na(df[[by]]))) stop('Dataframe contains NA values for "by"!')
  if(length(df[[by]]) != length(unique(df[[by]]))) stop('Dataframe contains duplicate observations of ', by, '!')
  
  ### Instead of raster::subs (which is pretty slow), just make a 
  ### vector of all the new values by joining the dataframe to a
  ### new dataframe with complete cell_id column, then replace all 
  ### raster values at once, very fast!
  out_rast <- get_mol_rast()
  df1 <- data.frame(cell_id = 1:ncell(out_rast)) %>%
    dt_join(df, by = 'cell_id', type = 'left')
  values(out_rast) <- df1[[which]]
  
  ### this keeps the layer name from the base rast... swap with "which"
  names(out_rast) <- which

  if(!is.null(xfm)) {
    if(class(xfm) != 'function') stop('xfm argument must be a function!')
    out_rast <- xfm(out_rast)
  }
  
  if(ocean_mask) {
    out_rast <- out_rast %>%
      terra::mask(terra::rast(here('_spatial/ocean_area_mol.tif')))
  }
  return(out_rast)
}

# comp_terra <- function(r1, r2, stopiffalse = TRUE) {
#   ### analogous to compareRaster(values = TRUE) for two
#   ### terra::SpatRaster objects
#   val_check <- values(r1) == values(r2)
#   val_check <- val_check[!is.na(val_check)]
#   
#   result <- all(val_check)
#   if(stopiffalse & !result) {
#     stop('Not all values are same between ', names(r1), ' and ', names(r2), '!!!')
#   }
#   return(result)
# }

#####################################################.
####   Helper functions for vulnerability data   ####
#####################################################.
get_spp_traits <- function() {
  ### this grabs all the traits by species, after downfilling and gapfilling
  trait_dir <- here_anx('../spp_vuln_framework/2_gapfill_traits')
  traits_gf_df <- data.table::fread(file.path(trait_dir, 'spp_up_down_gapfill_spp_traits.csv')) %>%
    full_join(data.table::fread(file.path(trait_dir, 'spp_up_down_gapfill_levels.csv'))) %>%
    full_join(data.table::fread(file.path(trait_dir, 'spp_up_down_gapfill_traits.csv'))) %>%
    full_join(data.table::fread(file.path(trait_dir, 'spp_up_down_gapfill_taxa.csv'))) %>%
    select(-starts_with('gf_'))
  return(traits_gf_df)
}

get_spp_vuln <- function() {
  ## to reassemble all scores incl components:
  vuln_dir <- here_anx('../spp_vuln_framework/3_vuln_score_traits')
  spp_vuln_scores_all <- data.table::fread(file.path(vuln_dir, 'spp_vuln_from_traits_all_scores.csv')) %>%
    left_join(data.table::fread(file.path(vuln_dir, 'spp_vuln_from_traits_str.csv')),
              by = 'vuln_str_id') %>%
    left_join(data.table::fread(file.path(vuln_dir, 'spp_vuln_from_traits_tx.csv')),
              by = 'vuln_tx_id') %>%
    select(-vuln_tx_id, -vuln_str_id) %>%
    select(species, taxon, stressor, score = vuln)
  spp_vuln_scores_fixed <- spp_vuln_scores_all %>%
    mutate(stressor = case_when(stressor == 'oa' ~ 'ocean_acidification',
                                stressor == 'uv' ~ 'uv_radiation',
                                stressor == 'slr' ~ 'sea_level_rise',
                                TRUE ~ stressor))
  return(spp_vuln_scores_fixed)
}

get_fe_traits <- function() {
  fe_trait_codes <- data.table::fread(here('_output/func_entities/fe_trait_codes.csv'))
  fe_traits <- data.table::fread(here('_output/func_entities/fe_species.csv')) %>%
    left_join(fe_trait_codes, by = c('fe_code', 'n_spp')) %>%
    select(-fe_code, -n_spp)
}

#####################################################################.
#### assemble overall species taxonomy from extracted WoRMS data ####
#####################################################################.
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


#####################################################.
####  Resolve names using WoRMS API fuzzy match  ####
#####################################################.

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

####################################################.
####          Raster support functions          ####
####################################################.

mask_ocean <- function(r) {
  ### for a raster in 10 km Mollweide, mask it to just
  ### ocean cells
  ocean_r <- raster(here('_spatial/ocean_area_mol.tif'))
  if(!compareRaster(r, ocean_r)) {
    warning('Cannot mask raster - CRS does not match')
    return(NULL)
  }
  x <- mask(r, ocean_r)
  return(x)
}

get_qtiles <- function(x) {
  q_vec <- c(0.5, 0.9, 0.95, 0.99, 0.999, 1.0)
  if(str_detect(tolower(class(x)), 'raster')) {
    x <- values(x)
  }
  quantile(x, q_vec, na.rm = TRUE)
}


####################################################.
####         Pooled variance functions          ####
####################################################.

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

### test that the function returns the correct variance value
# set.seed(42)
# n_vec <- sample(1:20, size = 123, replace = TRUE)
# x_list <- list()
# for(x_i in seq_along(n_vec)) { ### x_i <- 1
#   x_list[[x_i]] <- rnorm(mean = sample(5 * c(1:10), size = 1),
#                        sd   = sample(.5 * c(1:10), size = 1),
#                        n = n_vec[x_i])
# }
# 
# ### initialize values for first term
# mean_vec <- sapply(x_list, mean)
# sdev_vec <- sapply(x_list, sd)
# n_vec    <- sapply(x_list, length)
# 
# var_out <- iterated_pooled_var(mean_vec, sdev_vec, n_vec)
# 
# var_check <- var(unlist(x_list))
# 
# abs(var_check - var_out) < 1e-10

####################################################.
####     Functional vulnerability functions     ####
####################################################.

calc_fe <- function(n_spp) {
  k <- n_spp - 1
  fv <- 0.5^k
} 

calc_spp_cell_fe <- function(spp_cells, spp_fe) {
  ### parallelize this across smaller chunks to keep group_by from crashing 
  ### everything - but not for every cell individually!  
  ### Set up 1000 different cell groups across the 100k(ish) cells in the chunk
  ### to divide work between dplyr and parallel...
  cell_id_df <- data.frame(cell_id = spp_cells$cell_id %>% unique()) %>%
    mutate(cell_gp = rep(1:250, length.out = n()))
  cell_gps <- cell_id_df$cell_gp %>% unique()
  
  fe_df <- spp_cells %>%
    oharac::dt_join(spp_fe %>% select(species, fe_id), 
                    by = 'species', type = 'left')
  
  ### use the number of observations to limit the number of cores... more
  ### observations, fewer cores!
  n_cores <- ceiling(25 / ceiling((nrow(fe_df) / 8e7)))
  
  fv_list <- parallel::mclapply(cell_gps, mc.cores = n_cores,
                FUN = function(gp) { 
                  ### gp <- 5
                  cell_ids <- cell_id_df %>% filter(cell_gp == gp) %>% .$cell_id
                  
                  x <- fe_df %>%
                    filter(cell_id %in% cell_ids) %>%
                    data.table() %>% 
                    ### data.table syntax for: group_by/mutate (see below)
                    .[, ':='(n_spp = length(unique(species)),
                             n_fe  = length(unique(fe_id))),
                      by = .(cell_id)] %>%
                    .[, ':='(n_spp_fe = length(unique(species))),
                      by = .(cell_id, fe_id)] %>%
                    .[, ':='(fv = calc_fe(n_spp_fe)),
                      by = .(cell_id, fe_id)]
                  
                  return(x)
                })
  if(check_tryerror(fv_list)) {
    stop('Try-error found when calculating species/cell functional vulnerability!')
  }
  fv_df <- data.table::rbindlist(fv_list)
  return(fv_df)
}

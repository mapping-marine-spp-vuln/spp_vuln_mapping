# Organization of repository

All scripts and finalized outputs are contained within the Git-tracked repository.  Most external data and many intermediate files generated in the scripts are stored in an external directory on an NCEAS server.  To reproduce the entire analysis from scratch, you will need to create a similar external file structure (see below for details).

## Directory structure: Git-tracked

### Underscored directories: data and resources

* `_raw_data` includes all the data novel to this project, including expert-provided traits for a range of taxa, and sensitivity/adaptive capacity scores for traits across a range of stressors.
* `_data` includes data processed from various external sources (stressor distributions, species distributions, species vulnerability scores) used within the overall analysis.  This is generally *not* the original data.
* `_output` includes final products - rasters of impact by various weighting methods.

### Numbered directories: scripts

* `1_setup` includes scripts to process raw data
    * `functional_traits` directory includes scripts to collect data from FishBase/SealifeBase and the vulnerability framework project; assemble these data and impute missing values; and then assign species to functional entities based on these data.
    * `species` directory includes scripts to gather species range information and project to the base raster:
        * `iucn`: scrape information from IUCN API and then rasterize IUCN species range shapefiles, including depth cutoffs
        * `aquamaps`: clean/standardize AquaMaps species names, reproject half-degree raster to 10-km Mollweide, and clip to appropriate depth zones where necessary
        * For both IUCN and AquaMaps species, resulting presence maps are saved to an external location as .csvs with presence/cell ID (IUCN) or probability/cell ID (AquaMaps).
    * `1`stressors` directory includes scripts to reproject stressor layers from existing sources (e.g., Halpern et al. 2019 update, OHI, new raw data) to the 10 km Mollweide base raster.  
        * Most layers are "uniform exposure" stressors (i.e., all species exposed equally, though not necessarily equally vulnerable). For these, a single rescaled intensity map suffices.
        * Targeted fishing is species-specific, based on Watson 2018 data.  For these, a single intensity map will not be sufficient; instead, a map (in .csv form, per cell ID) of NPP-normalized fishing catch is generated, then these are normalized by either (a) the highest 90th percentile cell value seen across all stocks, or (b) the 99th percentile cell value for the individual stock, whichever is lower.
        * Thermal tolerance is species-specific, based on thermal envelopes from AquaMaps (or estimated from IUCN range in a similar manner to AquaMaps).  Similar to targeted fishing, a stressor intensity map is generated for each species and stored as a .csv (values per pixel by cell ID).
        * Bycatch is split into benthic and pelagic, and species are assigned to one or both bins based on water column position (benthic, pelagic, benthopelagic, reef).
* `2_process_spp_vuln_and_impacts` contains scripts to gather species rangemaps, functional entity memberships, vulnerability scores, and stressor maps to produce rasters of vulnerability and impact.  Maps of impacts by stressor or by taxon are saved in an external directory; aggregated cumulative impact maps are saved in the `_output` directory.
* `3_process_hab_vuln_and_impacts` contains scripts to calculate impacts according to a representative habitat method (e.g., Halpern et al 2008 or 2019).  As for the species-method impact maps, maps of impacts by stressor are saved in an external directory; aggregated cumulative impact maps are saved in the `_output` directory.
* `4_analysis` contains scripts to explore patterns in the results that are discussed in the manuscript.
* `5_ms_figs` contains scripts to generate the figures for the manuscript.

#### Running the scripts

Directories with an underscore prefix (e.g., `_data`) do not contain scripts, with one exception.  In `_spatial` there is a script, `set_up_master_rasters.Rmd`.  If you wish to reprocess this entire analysis in a different projection, resolution, etc., you may wish to investigate this script.  Otherwise, the resulting rasters from this script are already available in `_spatial` directory.

The numbered directories contain scripts (also numbered) that should be run in order to reproduce the results.

### Other directories

These are not critical to reproducing the results and can be generally ignored.

* `figs` includes auto-generated figures for other scripts; these are for data exploration purposes, not for inclusion in manuscript.
* `int` includes intermediate files saved during computationally intensive scripts to avoid re-processing and to save time.  These files may be used across multiple scripts.
* `explore` is old explorations and attempts at processing or analysis; this will probably be deleted.

## Directory structure: external (not Git-tracked)

We prefer to not host external datasets within our GitHub repository, and instead direct researchers to access the data from original sources (and with appropriate citations).  This ensures researchers are using the latest versions of data.  External files (external datasets and large intermediate analysis outputs) are hosted on Mazu server at NCEAS, in a directory `/home/shares/ohi/spp_vuln/spp_vuln_mapping/`.  To facilitate accessing these files, I've created a function analogous to `here::here()` that instead points to this "annex" directory:

```
here_anx <- function(f = '', ...) { 
  ### create file path to git-annex dir for project
  f <- paste(f, ..., sep = '/')
  f <- stringr::str_replace_all(f, '\\/+', '/')
  f_anx <- sprintf('/home/shares/ohi/spp_vuln/spp_vuln_mapping/%s', f)
  return(f_anx)
}
```

This function is available in `common_fxns.R` script in the project root directory.  Update the `f_anx` location to point to your preferred root location external to the Github repository.

I will not go into detail on the external directories for storing intermediate outputs; those can be found within the various scripts.  However, several scripts access data from external locations - you may wish to recreate a similar file structure or rewrite the script to accommodate a different file structure as you prefer.

### Species traits

Species traits (from Butt et al. 2022 Trait-Based Framework) are located in a parallel directory to the `spp_vuln_mapping` annex:

`/home/shares/ohi/spp_vuln/spp_vuln_framework/2_gapfill_traits`

In this directory are four files created in the workflow of the https://github.com/mapping-marine-spp-vuln/spp_vuln_framework project:

* `spp_up_down_gapfill_spp_traits.csv`: numeric codes for species, trait, and match level; also probability of this trait for this species
* `spp_up_down_gapfill_levels.csv`: match-code to match rank text
* `spp_up_down_gapfill_taxa.csv`: taxa-code to taxa/species text
* `spp_up_down_gapfill_traits.csv`: trait-code to category/trait/value text

See the `get_spp_traits()` function in `common_fxns.R`.

### Species vulnerability scores

Species vulnerability scores (from Butt et al. 2022 Trait-Based Framework) are located in a parallel directory to the `spp_vuln_mapping` annex:

`/home/shares/ohi/spp_vuln/spp_vuln_framework/2_gapfill_traits`

In this directory are four files created in the workflow of the https://github.com/mapping-marine-spp-vuln/spp_vuln_framework project:

* `spp_vuln_from_traits_all_scores.csv`: data frame of codes for taxon (by species) and stressor, along with sensitivity, general adaptive capacity, specific adaptive capacity, exposure, and combined vulnerability scores for each species (per taxon ID)
* `spp_vuln_from_traits_str.csv`: stressor code to stressor name
* `spp_vuln_from_traits_tx.csv`: taxonomic code to species/taxon name

See the `get_spp_vuln()` function in `common_fxns.R`.

### AquaMaps data

AquaMaps species distribution data (from aquamaps.org) are located in a parallel directory to the `spp_vuln_mapping` annex, as stored in the variable `am_dir` in `common_fxns.R`

`/home/shares/ohi/spp_vuln/aquamaps_2021`

In this directory are files sent from AquaMaps (direct communication with Kristin Kaschner).  Key files include (names may differ depending on how you receive them!):

* `hcaf_v7.csv`: half-degree cell authority file (HCAF) with environmental and geographic information on each half-degree cell
* `HCAFMetadata_v7.xls`: metadata for the HCAF info                         
* `ver10_2019_speciesoccursum_iucn.csv`: summary information on each species included in the dataset, including AquaMaps species ID, taxonomic info, occurrence cells (recommended to drop species with fewer than 10 occurrence cell records), and IUCN species ID.
* `ver10_2019_hspen.csv`: Environmental preference envelope data for species, used to determine distributions
* `species_envelope_summary.csv`: a summary of the environmental preference envelopes in tidy format      
* `hcaf_species_native.csv`: species ID to cell ID mapping, including probability of occurrence (relative environmental suitability)
* `hcaf_species_native_clean.csv`: a cleaner version of the species ID to cell ID mapping, using integer cell IDs instead of character, and dropping unneeded columns, for smaller file size and quicker read times (code to create this file is found in `1_setup/species/aquamaps/2_map_aquamaps_to_moll.Rmd`)

See helper functions in `common_fxns.R`:

* `get_am_spp_info()`
* `get_am_spp_envelopes()`
* `get_hcaf_info()`
* `get_hcaf_rast()`
* `map_to_hcaf()`
* `get_am_spp()`

### IUCN data

IUCN species range maps, downloaded from the [IUCN Red List Spatial Data Download site](https://www.iucnredlist.org/resources/spatial-data-download) (with a separate file for birds available upon request from BirdLife International, see note on previous line), are stored separately from the species vulnerability mapping project.  For our analysis, we stored them on Mazu server at NCEAS.  You can store them anywhere, but in script `1_setup/species/iucn/4_rasterize_spp_shps.Rmd` you will need to update these locations:

* dir_bli <- '/home/shares/ohi/git-annex/globalprep/_raw_data/birdlife_intl/d2021'
* dir_shp <- '/home/shares/ohi/git-annex/globalprep/_raw_data/iucn_spp/d2021-3'

All other data is retrieved from the [IUCN Red List API v3](https://apiv3.iucnredlist.org/api/v3/docs); much of this is stored externally as well, but you can find the external directory structure in the scripts, relative to the `here_anx()` location.

### Stressor data

See individual stressor processing scripts in `1_setup/stressors` for information on original sources of stressor data for most stressor types.  In some cases, stressor inputs are based on layers already processed for the Ocean Health Index or other NCEAS projects.  In these cases, please contact us (ohara@nceas.ucsb.edu or cohara@ucsb.edu) for the stressor data.

The resulting stressor maps, at 10 km Mollweide CRS, are all contained in the GitHub repository.




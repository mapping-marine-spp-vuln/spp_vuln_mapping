To Do:

- identify a few candidate stressors with data in Halpern et al. 2019
    - shipping (for ship strikes)
    - SST or OA
    - fishing? this seems like it needs a discussion
    - land-based pollutants
- calculate species sensitivity to the relevant stressors, both direct matches and gapfill
- identify species in IUCN data that match species in our vulnerability dataset; separate into direct and gapfilled at genus/family/order levels.

Calculate CHI on species range:
- calculate impact across range: sum of [species presence] $\times$ [stressor vulnerability] $\times$ [stressor intensity] across range
    - how to account for range sizes?  normalize by [range area]^_a_^ where _a_ is some exponent between 0 and 1?

Calculate map of impacted species:
- for each spp, as in BD/CHI, calculate some threshold of impact, i.e., [stressor vulnerability] $\times$ [stressor intensity], and map the impacted footprint.
- stack these to provide an "impact richness" map.  Normalize by number of spp present overall (again, like BD/CHI)
- for several areas with high and low proportional impact richness, take a patch (1+ cells) and calculate, for all spp, the cumulative impact from multiple stressors in that cell.  Plot this as stacked bars, one for each species present, with stacks representing contribution of multiple stressors on the species.  This gives a distribution of cumulative human impacts within the cell.

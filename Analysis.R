# =============================================================================
# analysis.R
# Consolidated longitudinal gut virome/microbiome analysis.
#
# This merges what used to be three separate scripts (master_tables.R,
# visualizations.R, taxonomy.R) into one file that runs top-to-bottom with
# no cross-file dependency ordering to get wrong. Source this ENTIRE file
# (RStudio: Source button, or Ctrl+Shift+Enter / Cmd+Shift+Return) after a
# clean Restart R - do not paste chunks into the console, and do not run
# only part of the file, since later PARTs depend on objects and helper
# functions (e.g. antibiotic_marker()) defined earlier in this same file.
#
#   PART 1: MASTER TABLES       - load + join + normalize + persistence
#   PART 2: VISUALIZATIONS      - diversity, abundance, persistence figures,
#           plus (2.6) the Aim A formal trend test (per-participant Kendall
#           trend test, richness/Shannon vs timepoint, BH-corrected)
#   PART 3: TAXONOMY            - viral/bacterial composition, CV trajectories
#           (3.6/3.7 most variable, 3.6B/3.7B least variable/most stable)
#   PART 4: FUNCTIONAL ANNOTATION - Phold gene function, lifestyle, AMGs,
#           plus (4.5) the Aim B per-function Fisher tests, temperate vs virulent (vOTU-level, BH)
#           and (4.6) the AMG pre/post-antibiotic-window Wilcoxon test
#   PART 5: CO-OCCURRENCE       - exploratory within-participant vOTU pairs
#   PART 6: VOLATILE vOTU FUNCTIONAL CONTEXT - what do the CV-selected vOTUs
#           from Part 3 actually carry, functionally (Part 4)?
#   PART 7: PHAGE:HOST RELATIONSHIP DYNAMICS - predicted-host pairing,
#           paired trajectories, correlation, antibiotic-window comparison
#   PART 8: COMMUNITY RESILIENCE - Aitchison-distance-from-baseline
#           trajectory, the concrete "does the community recover after
#           L1's antibiotic perturbation" metric Aim C was missing
#   PART 9: NETWORK STABILITY - L1's vOTU co-occurrence network, pre vs
#           post antibiotic: node/edge counts, density, mean edge
#           strength, and Jaccard overlap of the edge sets, plus (9.6) a
#           min_shared=3 sensitivity check and (9.7) a formal pre/post
#           diversity+abundance Wilcoxon test (viral and bacterial)
#   PART 10: DESCRIPTIVE COMMUNITY COMPOSITION - purely descriptive (no
#            tests): most-ABUNDANT (not CV-selected) genus/vOTU
#            trajectories, plus a top-20-taxa abundance heatmap per
#            participant/domain, at both genus (10.3) AND species (10.3B)
#            level for the bacterial side
#   PART 11: ADDITIONAL ANALYSES (NEW 2026-10-07) - PCoA ordination (all
#            participants + L1-only), within- vs between-person distance,
#            consecutive-sample membership turnover (viral vs bacterial),
#            richness-vs-sequencing-depth check, which bacterial genera
#            responded to L1's antibiotic (exploratory), temperate-share
#            pre/post antibiotic, and predicted host genus breakdown
#
# Each PART keeps its own internal section numbering (Part N -> N.1, N.2...)
# so it still matches earlier discussion/notes. Dependency chain (each part
# only ever ADDS columns/objects, never removes, so this order is safe):
#   Part 1 (raw tables; adds sample_dates, collection_date on viral_long/
#            bacterial_long, l1_antibiotic_marker_date)
#     -> Part 2 (adds antibiotic_marker()/antibiotic_marker_date()/
#                 integer_breaks() helpers, detection_freq_combined)
#       -> Part 3 (adds deepest_rank/taxon_label to viral_long, candidate_votus,
#                   candidate_genera, and the zero-filled genus_reads_full/
#                   species_reads_full/votu_reads_full tables Part 10 reuses)
#         -> Part 4 (adds phold_votu; needs viral_long's TYPE/lifestyle column
#                     from Part 1 and antibiotic_marker() from Part 2)
#           -> Part 5 (adds clr_transform_matrix(), cooccurrence_filtered;
#                       needs viral_long/abundance_wide from Part 1)
#             -> Part 6 (needs candidate_votus + taxonomy from Part 3, and
#                         phold_votu from Part 4)
#             -> Part 7 (needs viral_long/bacterial_long from Part 1,
#                         antibiotic_marker() from Part 2, and reuses Part 5's
#                         framing/caveats around small-n exploratory stats;
#                         adds bacterial_long$bracken_rel_abundance)
#             -> Part 8 (needs viral_long/bacterial_long/abundance_wide from
#                         Part 1, antibiotic_marker() from Part 2, and
#                         clr_transform_matrix() from Part 5)
#             -> Part 9 (needs viral_long/abundance_wide from Part 1, and
#                         clr_transform_matrix() from Part 5; 9.7 also needs
#                         viral_diversity/bacterial_diversity/
#                         viral_total_abundance from Part 2)
#               -> Part 10 (needs genus_reads_full/species_reads_full/
#                            votu_reads_full from Part 3, antibiotic_marker_date()
#                            from Part 2)
#                 -> Part 11 (needs viral_long/bacterial_long/abundance_wide
#                              from Part 1, genus_reads_full from Part 3,
#                              clr_transform_matrix() from Part 5,
#                              get_participant_clr() from Part 8,
#                              window_wilcoxon_test() from Part 9,
#                              antibiotic_marker_date()/
#                              participant_date_scales()/integer_breaks() from
#                              Part 2)
# So PART 1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 7 -> 8 -> 9 -> 10 is the only valid order.
# =============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)
library(forcats)
library(ggplot2)

# Console tables: never truncate columns or long strings (several tables were
# cut off at 80 characters in earlier runs).
options(tibble.width = Inf, pillar.width = Inf, width = 250)

# NEW 2026-09-29: ggh4x, for facetted_pos_scales() - gives each participant's
# panel its own exact x-axis date breaks in the handful of figures that show
# all 3 participants side by side (see Part 2 for why a single shared set of
# date breaks isn't correct here). Installed automatically below if missing,
# so a fresh machine doesn't need a manual install.packages() step first.
if (!requireNamespace("ggh4x", quietly = TRUE)) {
  message("Installing 'ggh4x' (one-time - not previously used in this project)...")
  install.packages("ggh4x")
}
library(ggh4x)

# -----------------------------------------------------------------------
# COLOR PALETTE (NEW 2026-09-30) - one small set of colorblind-safe hex
# values, reused everywhere by MEANING rather than picked per-plot. This
# replaces what had become ~10 unrelated named colors (steelblue,
# chocolate, seagreen, darkslategray, midnightblue, darkorange3, ...)
# scattered across the script with no consistent logic, which made
# domain (viral vs bacterial) hard to track visually from one figure to
# the next and was never checked for colorblind-safety. Values are a
# validated colorblind-safe categorical palette (fixed hue order - the
# order itself is the CVD-safety mechanism, so it's never reshuffled
# per-plot; each adjacent pair stays distinguishable under simulated
# protanopia/deuteranopia). The two domain colors also double as the
# base hue for each domain's sequential heatmap ramp (viral -> blue
# "viridis", bacterial -> orange "magma", both already perceptually
# uniform, colorblind-safe continuous scales), so a reader associates
# the same hue family with the same domain everywhere in the thesis,
# not just in the line charts.
# -----------------------------------------------------------------------
COL_VIRAL     <- "#2a78d6"  # blue    - every viral-domain series/marker
COL_BACTERIAL <- "#eb6834"  # orange  - every bacterial-domain series/marker
COL_RICHNESS  <- "#1baf7a"  # aqua    - viral richness (kept distinct from viral diversity/abundance on shared figures)
COL_AMG       <- "#eda100"  # yellow  - AMG-carrying vOTU trend
COL_ABUNDANCE <- "#e87ba4"  # magenta - viral total-abundance trend
COL_MARKER    <- "#e34948"  # red     - antibiotic-window marker line (an event, not a data series)
COL_DIV_LOW   <- "#e34948"  # red     - diverging scale, negative pole
COL_DIV_MID   <- "#f0efec"  # neutral gray - diverging scale midpoint (0 correlation)
COL_DIV_HIGH  <- "#2a78d6"  # blue    - diverging scale, positive pole
# Within-domain ORDINAL shading for "most variable" vs "least variable/
# stable" trend panels - same hue family per domain, dark vs light step,
# since "most/least variable" is a position on one spectrum (an ordinal
# distinction), not a separate categorical identity.
COL_VIRAL_DARK      <- "#184f95"  # most variable vOTUs
COL_VIRAL_LIGHT     <- "#5598e7"  # least variable (most stable) vOTUs
COL_BACTERIAL_DARK  <- "#b34d1f"  # most variable genera/species
COL_BACTERIAL_LIGHT <- "#f0a878"  # least variable (most stable) genera/species
# Fixed per-participant identity colors - used anywhere D3/G2/L1 are
# distinguished by fill/colour within the SAME panel (rather than split
# across facets, which needs no color coding at all).
PARTICIPANT_COLORS <- c(D3 = "#1baf7a", G2 = "#e87ba4", L1 = "#4a3aa7")


# #############################################################################
# PART 1: MASTER TABLES
# (originally master_tables.R)
# Build master long-format tables (viral, bacterial) for the longitudinal
# gut virome/microbiome thesis analysis.
# #############################################################################

# -----------------------------------------------------------------------
# 1.0. Paths - local Windows copy (transferred via scp from Setonix)
# -----------------------------------------------------------------------
BASE <- "C:/Users/love0371/OneDrive/Documents/Flinders_Phage_Gut_Microbiome_Project/Bioinformatics"

# Where the two large co-occurrence heatmaps (Parts 5 and 9) get exported
# to at full size (NEW 2026-09-29) - RStudio's Plots-pane preview is too
# small for a heatmap with dozens of vOTU labels on each axis to ever
# read cleanly, however small the font is set; these are saved as actual
# image files sized to the number of labels instead. Created on first run
# if it doesn't already exist.
FIGURES_DIR <- file.path(BASE, "results", "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------
# 1.1. VIRAL MASTER TABLE
# -----------------------------------------------------------------------

## 1.1a. vOTU abundance (wide: contig x sample) -> long
abundance_wide <- read_tsv(
  file.path(BASE, "results/abundance/votu_abundance_table.tsv"),
  show_col_types = FALSE
)

abundance_long <- abundance_wide |>
  pivot_longer(
    cols      = -contig,
    names_to  = "sample_id",
    values_to = "abundance"
  ) |>
  separate(sample_id, into = c("participant", "timepoint"),
           sep = "\\.", remove = FALSE) |>
  mutate(timepoint = as.numeric(timepoint))

## 1.1b. Confidence tiers (one row per vOTU representative contig)
confidence_tiers <- read_tsv(
  file.path(BASE, "results/votu/confidence_tiers.tsv"),
  show_col_types = FALSE
)

## 1.1c. PhaBOX2 annotations - host, taxonomy, lifestyle (all keyed by Accession)
host <- read_tsv(
  file.path(BASE, "results/annotation/phabox/host_prediction.tsv"),
  show_col_types = FALSE
) |>
  rename(contig = Accession)

taxonomy <- read_tsv(
  file.path(BASE, "results/annotation/phabox/taxonomy.tsv"),
  show_col_types = FALSE
) |>
  rename(contig = Accession)

lifestyle <- read_tsv(
  file.path(BASE, "results/annotation/phabox/lifestyle.tsv"),
  show_col_types = FALSE
) |>
  rename(contig = Accession)

## 1.1d. Join everything into one long viral master table
viral_long <- abundance_long |>
  left_join(
    confidence_tiers |> rename(contig = representative_contig),
    by = "contig"
  ) |>
  left_join(host,      by = "contig") |>
  left_join(taxonomy,  by = "contig") |>
  left_join(lifestyle, by = "contig")

## 1.1e. Join diagnostics - run this BEFORE trusting the table.
## A high unmatched percentage on any of these means the join key
## (contig ID format) doesn't line up across files.
check_join <- function(df, col, label) {
  n_total   <- nrow(df)
  n_missing <- sum(is.na(df[[col]]))
  pct <- round(100 * n_missing / n_total, 1)
  message(sprintf("%-25s: %5d / %5d rows unmatched (%s%%)",
                  label, n_missing, n_total, pct))
}
message("=== Join diagnostics (viral_long) ===")
check_join(viral_long, "tier",  "confidence_tiers join")
check_join(viral_long, "Host",  "host_prediction join")
check_join(viral_long, "Genus", "taxonomy join")
check_join(viral_long, "TYPE",  "lifestyle join")
# If confidence_tiers join looks broken: inspect
# confidence_tiers$representative_contig directly - a value seen during
# development ("L1" instead of a full contig ID like "L1.9__contig_1313")
# looked truncated/wrong and may indicate a bug in
# scripts/assign_confidence_tiers.py rather than a problem with this join.
print(table(nchar(confidence_tiers$representative_contig) > 5))

# -----------------------------------------------------------------------
# 1.1f. CLEAN vOTU DISPLAY IDs (NEW 2026-09-29) - for every figure from
# here on, never for joins.
# -----------------------------------------------------------------------
# `contig` (== representative_contig from CD-HIT dereplication) is named
# "{participant}.{sample}__contig_{n}", where {participant} is whichever
# SAMPLE's assembly happened to contribute the representative sequence
# for that vOTU cluster - NOT necessarily the participant a given figure
# is showing that vOTU under. A vOTU detected in L1 can easily carry a
# "D3." prefix if D3's assembly produced the best representative
# sequence during dereplication, which reads as a labelling error to
# anyone looking at "L1's" figure, even though it isn't one. Fixed by
# giving every vOTU a clean, sequential, participant-agnostic display
# label (vOTU_001, vOTU_002, ...) used on every figure from here on. The
# real contig ID is never discarded - `votu_display_lookup` keeps the
# mapping back to it, which is worth including as a supplementary table
# in the thesis so a reader can still trace a specific vOTU_### back to
# its actual contig/assembly if needed.
votu_display_lookup <- viral_long |>
  distinct(contig, votu_id) |>
  arrange(contig) |>
  mutate(votu_display = sprintf("vOTU_%03d", row_number()))

message(sprintf("=== vOTU display-ID lookup built: %d vOTUs relabelled vOTU_001-vOTU_%03d for every figure (see votu_display_lookup for the real contig IDs behind each one) ===",
                nrow(votu_display_lookup), nrow(votu_display_lookup)))

viral_long <- viral_long |>
  left_join(votu_display_lookup, by = c("contig", "votu_id"))

# Canonical list of participants. Defined once, this early, because
# several later Parts (4, 5, 8, 9) each iterate "over every participant"
# and previously only Part 5 (Section 5.3) created this object - a script
# run top-to-bottom from a clean R session would hit Part 4 Section 4.6's
# AMG Wilcoxon test first and fail with "object 'all_participants' not
# found" (found 2026-09-29 during a full re-check; it likely only ever
# "worked" before because a stale copy was left over in memory from an
# earlier partial/non-fresh run, not because the script actually ran
# clean top-to-bottom).
all_participants <- unique(viral_long$participant)

# -----------------------------------------------------------------------
# 1.2. BACTERIAL MASTER TABLE (Bracken)
# -----------------------------------------------------------------------

bracken_files <- list.files(
  file.path(BASE, "results/bacterial/bracken"),
  pattern    = "_bracken\\.tsv$",
  full.names = TRUE
)

read_one_bracken <- function(path) {
  sample_id <- str_remove(basename(path), "_bracken\\.tsv$")
  read_tsv(path, show_col_types = FALSE) |>
    mutate(sample_id = sample_id)
}

bacterial_long <- map_dfr(bracken_files, read_one_bracken) |>
  separate(sample_id, into = c("participant", "timepoint"),
           sep = "\\.", remove = FALSE) |>
  mutate(timepoint = as.numeric(timepoint)) |>
  rename(bacterial_taxon = name)

## 1.2a. Exclude phage/viral contamination (found 2026-09-25). Kraken2's
## standard database includes viral genomes alongside bacterial ones, so
## Bracken's "bacterial" species-level output isn't purely bacterial.
## Confirmed 28/573 distinct taxa (0.553% of total reads) are actually
## bacteriophage taxa - caught via ICTV phage genus naming convention
## (first word ending in "virus", e.g. "Endlipuvirus intestinihominis" -
## two of these genera, Endlipuvirus and Carjivirus, are also literally
## in the vOTU catalogue's own PhaBOX2 Genus column) or "phage" appearing
## in the species name (e.g. "Parabacteroides phage YZ-2015b"). This
## MUST happen here, before bacterial_long is used anywhere downstream
## (Part 1 Section 1.4's persistence comparison, Part 3's genus
## composition), not as a later filter - a narrower version of this
## exclusion used to live only in taxonomy.R (matching just
## "^uncultured phage"), which caught 1 of the 28 contaminant taxa and
## left the rest silently counted as "bacterial" in the
## persistence/composition results.
bacterial_long <- bacterial_long |>
  filter(
    !str_detect(word(bacterial_taxon, 1), "virus$"),
    !str_detect(bacterial_taxon, "[Pp]hage")
  )

message(sprintf("Bacterial long table: %d rows from %d samples",
                nrow(bacterial_long), length(bracken_files)))

## 1.2b. Exclude non-bacterial (host/human) contamination (found 2026-09-29,
## flagged by the user noticing "Homo" among the least-variable bacterial
## genera). Kraken2's standard database bundles the human reference genome
## (GRCh38) alongside bacterial genomes by default, so a small fraction of
## "bacterial" Bracken calls can actually be human host DNA that wasn't
## filtered out at the read-classification stage. Same
## must-happen-before-any-downstream-use logic as the phage exclusion above
## (Section 1.2a) - this runs before Part 1 Section 1.4's persistence
## comparison and Part 3's genus composition, not as a later filter.
##
## Bracken/Kraken2 doesn't report a domain/kingdom column in this output
## (species-level rows only), so a systematic domain-level filter isn't
## possible without a new dependency (e.g. a local NCBI taxonomy lookup) -
## this instead excludes known non-bacterial genera by name, the same
## targeted approach already used for phage. Currently just Homo (the
## concrete case found); if the diagnostic below turns up other eukaryotic
## contaminants (e.g. other host/lab species), add their genus name here.
known_non_bacterial_genera <- c("Homo")

n_before_host_filter <- nrow(bacterial_long)
bacterial_long <- bacterial_long |>
  filter(!word(bacterial_taxon, 1) %in% known_non_bacterial_genera)
n_removed_host <- n_before_host_filter - nrow(bacterial_long)

message(sprintf(
  "Bacterial long table: removed %d rows (genera: %s) as non-bacterial (host) contamination",
  n_removed_host, paste(known_non_bacterial_genera, collapse = ", ")
))

# Diagnostic only, not an automatic filter: flags any remaining genus whose
# first word doesn't look like an ordinary binomial-style bacterial genus
# (capitalized Latin word) - a quick eyeball check so anything else unusual
# doesn't silently slip through the way Homo did. Review this list in the
# console output and add confirmed non-bacterial genera to
# known_non_bacterial_genera above.
message("=== Bacterial genera worth a visual double-check (unusual first-word pattern) ===")
bacterial_long |>
  distinct(bacterial_taxon) |>
  mutate(first_word = word(bacterial_taxon, 1)) |>
  filter(!str_detect(first_word, "^[A-Z][a-z]+$")) |>
  distinct(first_word) |>
  print(n = Inf)

# -----------------------------------------------------------------------
# 1.2c. SAMPLE COLLECTION DATES
# -----------------------------------------------------------------------
# Real calendar collection dates (supplied 2026-09-29), used from here on
# as the x-axis for any plot tracking change over time. `timepoint`
# (sample order: 1, 2, 3...) compresses/stretches real-world gaps between
# visits into equal-width steps, which is especially misleading for L1
# (gaps of 14-49 days between samples) right around the antibiotic window.
# `timepoint` is KEPT and still used everywhere for windowing/ordering
# logic (pre/post-antibiotic split, baseline definitions, CV floors, etc.)
# since those are inherently about sample ORDER relative to the antibiotic
# course, not the calendar - only plotting x-axes change to
# collection_date, and only where a genuine time-series is being shown
# (see the per-plot notes at each PART for which plots did/didn't move to
# date axes and why).

sample_dates <- tibble(
  sample_id = c(
    paste0("L1.", 1:9),
    paste0("G2.", 1:8),
    paste0("D3.", 1:10)
  ),
  collection_date = as.Date(c(
    "2025-04-03", "2025-04-17", "2025-05-02", "2025-05-15", "2025-06-17",
    "2025-07-11", "2025-08-06", "2025-09-26", "2025-10-26",
    "2025-10-12", "2025-10-26", "2025-11-10", "2025-11-23", "2025-12-07",
    "2026-01-01", "2026-01-25", "2026-02-09",
    "2026-01-04", "2026-01-20", "2026-02-01", "2026-02-15", "2026-03-01",
    "2026-03-15", "2026-03-29", "2026-04-12", "2026-04-28", "2026-05-10"
  ))
)

viral_long     <- viral_long     |> left_join(sample_dates, by = "sample_id")
bacterial_long <- bacterial_long |> left_join(sample_dates, by = "sample_id")

message(sprintf(
  "Sample dates joined - unmatched sample_id rows (should both be 0): viral_long = %d, bacterial_long = %d",
  sum(is.na(viral_long$collection_date)), sum(is.na(bacterial_long$collection_date))
))

# L1's antibiotic course fell between samples L1.4 and L1.5 (used
# throughout as the numeric xintercept = 4.5 marker on timepoint-axis
# plots) - this is the same event's position on a real calendar axis,
# the midpoint between those two samples' actual collection dates, for
# the plots below that move to collection_date.
l1_antibiotic_marker_date <- mean(
  sample_dates$collection_date[sample_dates$sample_id %in% c("L1.4", "L1.5")]
)

# -----------------------------------------------------------------------
# 1.3. NORMALIZATION - GMPR (Geometric Mean of Pairwise Ratios)
# -----------------------------------------------------------------------
# NOTE: votu_abundance_table.tsv holds CoverM *trimmed-mean coverage*
# (continuous - see trim_min/trim_max in the coverm rule), NOT raw read
# counts.
#
# NORMALIZATION METHOD: this matrix is very sparse - some samples have
# only a handful of non-zero vOTUs out of the 151-vOTU catalogue (checked empirically). Standard
# median-of-ratios (DESeq2-style, using one fixed cross-sample reference)
# is unstable here, since a sample's factor can end up based on a
# handful of ratios. GMPR (Chen et al. 2018) is designed for exactly
# this situation: it computes ratios pairwise, using only the features
# two samples actually share, rather than requiring a shared reference
# across all samples.

gmpr_size_factors <- function(count_matrix, min_shared = 3) {
  # count_matrix: contigs (rows) x samples (columns), numeric, >= 0
  n <- ncol(count_matrix)
  samples <- colnames(count_matrix)
  gmpr <- numeric(n)
  n_pairs_used <- integer(n)
  names(gmpr) <- samples
  names(n_pairs_used) <- samples
  
  for (j in seq_len(n)) {
    pairwise_ratios <- c()
    for (k in seq_len(n)) {
      if (k == j) next
      shared <- count_matrix[, j] > 0 & count_matrix[, k] > 0
      if (sum(shared) >= min_shared) {
        r <- median(count_matrix[shared, j] / count_matrix[shared, k])
        pairwise_ratios <- c(pairwise_ratios, r)
      }
    }
    n_pairs_used[j] <- length(pairwise_ratios)
    if (length(pairwise_ratios) == 0) {
      warning(sprintf(
        "Sample %s has no valid pairwise comparisons (min_shared=%d) - size factor set to NA",
        samples[j], min_shared
      ))
      gmpr[j] <- NA
    } else {
      gmpr[j] <- exp(mean(log(pairwise_ratios)))
    }
  }
  
  # Normalize so the geometric mean of all factors is 1 (interpretability)
  valid <- !is.na(gmpr) & gmpr > 0
  gmpr[valid] <- gmpr[valid] / exp(mean(log(gmpr[valid])))
  
  attr(gmpr, "n_pairs_used") <- n_pairs_used
  gmpr
}

abundance_matrix <- abundance_wide |>
  column_to_rownames("contig") |>
  as.matrix()

size_factors <- gmpr_size_factors(abundance_matrix)

message("=== GMPR: pairwise comparisons used per sample ===")
print(attr(size_factors, "n_pairs_used"))
message("Samples with few pairs (<5) may still have noisy size factors -")
message("worth treating their model results with extra caution.")

size_factor_df <- tibble(
  sample_id   = names(size_factors),
  size_factor = size_factors
)

viral_long <- viral_long |>
  left_join(size_factor_df, by = "sample_id")

# -----------------------------------------------------------------------
# 1.4. PERSISTENCE - READ-MAPPING-BASED (the correct persistence metric)
# -----------------------------------------------------------------------
# IMPORTANT: confidence_tiers$n_samples (from CD-HIT clustering) answers a
# different, much stricter question than "is this organism present in
# multiple samples" - it only counts a vOTU as multi-sample if the SAME
# representative contig independently assembled and clustered across >=2
# samples. That is an assembly/clustering-concordance metric, not a
# persistence metric, and treating it as persistence is what produced an
# earlier, biologically implausible "near-total viral turnover" finding.
#
# Read-mapping-based detection (CoverM abundance > 0 for viral, Bracken
# new_est_reads > 0 for bacterial) answers the biologically meaningful
# question instead, and is comparable across viral and bacterial data on
# the same basis. USE THIS SECTION'S OUTPUT for any persistence/stability
# claim in the thesis.

participant_n_samples <- viral_long |>
  distinct(participant, sample_id) |>
  count(participant, name = "total_samples")

## 1.4a. Viral: detected per contig per participant
mapping_based_persistence <- viral_long |>
  filter(abundance > 0) |>
  group_by(contig, participant) |>
  summarise(n_samples_detected = n_distinct(sample_id), .groups = "drop")

## 1.4b. Bacterial: detected per taxon per participant
bacterial_mapping_persistence <- bacterial_long |>
  filter(new_est_reads > 0) |>
  group_by(bacterial_taxon, participant) |>
  summarise(n_samples_detected = n_distinct(sample_id), .groups = "drop")

## 1.4c. Convert to detection FREQUENCY (0-1), since D3/G2/L1 have different
## total sample counts (10/8/9) - a raw n_samples_detected count isn't
## directly comparable across participants or between domains without
## this correction.
viral_freq <- mapping_based_persistence |>
  left_join(participant_n_samples, by = "participant") |>
  mutate(detection_freq = n_samples_detected / total_samples)

bacterial_freq <- bacterial_mapping_persistence |>
  left_join(participant_n_samples, by = "participant") |>
  mutate(detection_freq = n_samples_detected / total_samples)

detection_freq_combined <- bind_rows(
  viral_freq     |> transmute(participant, detection_freq, domain = "Viral"),
  bacterial_freq |> transmute(participant, detection_freq, domain = "Bacterial")
)

## 1.4d. Core (>=80% of samples) / intermediate / transient (<=20%) split -
## a concrete effect-size number to sit alongside the KS test in Part 2
## (a p-value alone says "different," not "how").
message("=== Detection frequency: core / intermediate / transient split ===")
detection_freq_combined |>
  mutate(category = case_when(
    detection_freq <= 0.2 ~ "transient (<=20% of samples)",
    detection_freq >= 0.8 ~ "core (>=80% of samples)",
    TRUE ~ "intermediate"
  )) |>
  count(domain, category) |>
  group_by(domain) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  arrange(domain, category) |>
  print(n = Inf)

# #############################################################################
# PART 2: VISUALIZATIONS
# (originally visualizations.R)
# Temporal visualisation of viral and bacterial community structure.
# Depends on Part 1's viral_long, bacterial_long, confidence_tiers,
# detection_freq_combined.
# #############################################################################

# -----------------------------------------------------------------------
# CAVEAT: `timepoint` is sample ORDER (1, 2, 3...), not actual elapsed
# time. L1's antibiotic course fell between timepoints 4 and 5, but if
# real collection dates are available, swap them in here - the ordinal
# axis can visually compress/stretch real gaps between samples,
# especially relevant for L1's irregular intervals.
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# Helper: mark L1's antibiotic window (between timepoints 4 and 5).
# Part 3 (taxonomy) also calls this helper, which is why Part 2 must run
# before Part 3.
# -----------------------------------------------------------------------
antibiotic_marker <- function(participant_filter = "L1", x_position = 4.5) {
  list(
    geom_vline(
      data = data.frame(participant = participant_filter, xintercept = x_position),
      aes(xintercept = xintercept),
      linetype = "dashed", colour = "firebrick", linewidth = 0.6
    )
  )
}

# Date-axis counterpart to antibiotic_marker(), for plots that now use real
# collection dates (Part 1, Section 1.2c) instead of sample order - draws
# the same dashed line, just at l1_antibiotic_marker_date (the midpoint
# between L1.4's and L1.5's real collection dates) instead of the fixed
# numeric position 4.5.
antibiotic_marker_date <- function(participant_filter = "L1") {
  list(
    geom_vline(
      data = data.frame(participant = participant_filter, xintercept = l1_antibiotic_marker_date),
      aes(xintercept = xintercept),
      linetype = "dashed", colour = "firebrick", linewidth = 0.6
    )
  )
}

# Integer-only x-axis breaks, for plots that still use `timepoint` (sample
# order) as the x-axis rather than collection_date - e.g. the stacked-
# proportion composition plots, where uneven real-world date spacing would
# make bar widths/gaps misleading rather than informative. Without this,
# ggplot's default breaks algorithm can place fractional ticks (e.g. 2.5,
# 7.5) on a short integer sequence, which doesn't correspond to a real
# sample. Works per-panel under facet_wrap(scales = "free_x"), since
# ggplot2 calls a function passed to `breaks` separately for each panel's
# own data range.
integer_breaks <- function(x) {
  rng <- range(x, na.rm = TRUE)
  seq(ceiling(rng[1]), floor(rng[2]))
}

# -----------------------------------------------------------------------
# Helper: exact per-participant date breaks for the handful of figures
# that show all 3 participants side-by-side in one faceted image
# (facet_wrap(scales = "free_x")). NEW 2026-09-29, fixes "L1's collection
# dates are partially missing on the x-axis".
#
# What was happening: scale_x_date() with no explicit `breaks` falls
# back to an automatic "pretty" algorithm that picks evenly-spaced
# calendar intervals - it was never actually using the real sample dates
# as its tick positions. L1's samples are by far the most irregularly
# spaced (14-49 days between visits, vs D3/G2's closer-to-fortnightly
# rhythm), so "pretty" intervals landed on some of L1's real dates by
# chance and skipped others entirely - that's the "partially missing"
# pattern.
#
# The fix isn't just "pass every real date as a break", though: L1, G2
# and D3's actual collection-date windows overlap at the edges (L1's
# last sample, 2025-10-26, is the exact same date as G2's 2nd sample;
# G2 and D3 both have samples in Jan 2026), so one shared breaks vector
# would put a stray tick/date label - with no real sample underneath it
# - onto a neighbouring participant's panel wherever the ranges cross.
# ggh4x::facetted_pos_scales() gives each facet panel a genuinely
# independent x scale, so every panel shows exactly (and only) that
# participant's own real sample dates, however irregular the spacing.
participant_date_scales <- function() {
  by_participant <- split(sample_dates$collection_date,
                          sub("\\..*", "", sample_dates$sample_id))
  ggh4x::facetted_pos_scales(
    x = lapply(by_participant, function(d) {
      scale_x_date(breaks = sort(unique(d)), date_labels = "%d %b %y")
    })
  )
}

# -----------------------------------------------------------------------
# Helper: move every figure's title below the plot as a caption, per
# scientific-figure convention (title/legend text sits under the image,
# not above it the way ggplot2 defaults to). NEW 2026-09-29.
#
# Subtitles carried extra methodological detail (exactly which
# test/transform was used, sample-size caveats, etc.) that's useful
# context but too cluttered to leave sitting on the figure image itself
# - style_figure() drops the subtitle from the image entirely, but never
# throws the wording away: both the title and subtitle text are stashed
# in `figure_captions` (every figure appends one entry, in the order it
# was drawn; printed as a single list at the very end of this script) so
# the exact wording is on hand to paste into the real thesis figure
# legends without having to reconstruct it from the code.
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Helper: shrink + reformat correlation-heatmap axis text as the number
# of distinct vOTUs grows, so a heatmap with many categories (a
# participant with lots of qualifying vOTUs - L1 especially, given it
# has the most timepoints) stays as legible as possible on-screen. Used
# by the Part 5 and Part 9 co-occurrence heatmaps, both of which are ALSO
# exported to FIGURES_DIR at a size that scales with n_labels (passed as
# style_figure()'s own width/height args - see each call site) - even a
# small, non-overlapping font still needs enough physical canvas
# width/height to fit 40-60+ labels along an axis, which the RStudio
# Plots-pane preview can't give it. NEW 2026-09-29, fixes "the contig
# labels on the side are all overlapping".
# -----------------------------------------------------------------------
heatmap_axis_theme <- function(n_labels) {
  fs <- pmax(3, pmin(6, 260 / n_labels))
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = fs),
        axis.text.y = element_text(size = fs))
}

# -----------------------------------------------------------------------
# Shared base theme (NEW 2026-09-30) - every figure in this script now
# builds on this instead of a bare theme_thesis(), so font sizes, strip
# (facet label) styling, and legend placement read as one consistent
# system across ~55+ figures rather than whatever each plot's own
# theme() call happened to set. Still a thin wrapper around
# theme_thesis() - individual plots' own theme() calls (axis text
# angle, per-heatmap font scaling, italic strip text for taxon-name
# facets, etc.) are unaffected, since those are added on top of this via
# `+` and simply override the relevant element.
# -----------------------------------------------------------------------
theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      strip.text       = element_text(face = "bold", size = rel(0.95)),
      plot.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

figure_captions <- list()

# -----------------------------------------------------------------------
# style_figure() now also handles export (NEW 2026-09-30, folded into
# the same step rather than kept separate): every figure that passes
# through here is saved to FIGURES_DIR as a 300dpi PNG with a consistent
# naming scheme (fig01_..., fig02_..., in draw order, slugged from its
# own title), as well as being print()'d to the Plots pane/PDF device -
# previously only the two co-occurrence heatmaps got this treatment,
# which meant every other figure only ever existed at RStudio's small
# preview resolution unless exported by hand. This is what makes every
# figure in the thesis "publication size" by construction rather than a
# manual per-figure export step. width/height/dpi default to a standard
# single-panel figure size but can be overridden per call - the two
# co-occurrence heatmaps do this, since their legible size depends on
# how many vOTU labels are on each axis (see Parts 5 and 9).
# -----------------------------------------------------------------------
style_figure <- function(p, width = 8, height = 6, dpi = 300) {
  ttl <- p$labels$title
  sub <- p$labels$subtitle
  
  figure_captions[[length(figure_captions) + 1]] <<- list(title = ttl, subtitle = sub)
  
  fig_n <- length(figure_captions)
  slug <- if (is.null(ttl)) {
    "untitled"
  } else {
    ttl |>
      tolower() |>
      str_replace_all("[^a-z0-9]+", "_") |>
      str_replace_all("^_+|_+$", "")
  }
  
  styled <- p +
    labs(title = NULL, subtitle = NULL, caption = ttl) +
    theme(
      plot.caption          = element_text(hjust = 0, size = 10, margin = margin(t = 8)),
      plot.caption.position = "plot"
    )
  
  ggsave(
    file.path(FIGURES_DIR, sprintf("fig%02d_%s.png", fig_n, slug)),
    plot = styled, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE
  )
  
  styled
}

# -----------------------------------------------------------------------
# 2.1. VIRAL DIVERSITY OVER TIME (richness + Shannon)
# -----------------------------------------------------------------------
# Diversity indices use within-sample relative abundance, so size-factor
# normalization is not needed here (a constant per-sample scaling factor
# cancels out when computing proportions within that sample).

viral_diversity <- viral_long |>
  group_by(sample_id, participant, timepoint) |>
  summarise(
    richness = sum(abundance > 0),
    shannon  = {
      p <- abundance[abundance > 0] / sum(abundance)
      -sum(p * log(p))
    },
    .groups = "drop"
  ) |>
  left_join(sample_dates, by = "sample_id")

p_viral_shannon <- ggplot(viral_diversity, aes(x = collection_date, y = shannon)) +
  geom_line(colour = COL_VIRAL) +
  geom_point(size = 2, colour = COL_VIRAL) +
  antibiotic_marker_date() +
  facet_wrap(~ participant, scales = "free_x") +
  participant_date_scales() +
  labs(title = "Viral (vOTU) Shannon diversity over time",
       x = "Collection date", y = "Shannon diversity") +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_viral_richness <- ggplot(viral_diversity, aes(x = collection_date, y = richness)) +
  geom_line(colour = COL_RICHNESS) +
  geom_point(size = 2, colour = COL_RICHNESS) +
  antibiotic_marker_date() +
  facet_wrap(~ participant, scales = "free_x") +
  participant_date_scales() +
  labs(title = "Viral (vOTU) richness over time",
       x = "Collection date", y = "Number of vOTUs detected") +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(style_figure(p_viral_shannon))
print(style_figure(p_viral_richness))

# -----------------------------------------------------------------------
# 2.2. TOTAL NORMALIZED VIRAL ABUNDANCE OVER TIME
# -----------------------------------------------------------------------
# Here we DO need size-factor normalization, since we're comparing total
# magnitude across samples rather than within-sample proportions.

viral_total_abundance <- viral_long |>
  group_by(sample_id, participant, timepoint, size_factor) |>
  summarise(total_abundance = sum(abundance), .groups = "drop") |>
  mutate(normalized_abundance = total_abundance / size_factor) |>
  left_join(sample_dates, by = "sample_id")

p_viral_abundance <- ggplot(viral_total_abundance, aes(x = collection_date, y = normalized_abundance)) +
  geom_line(colour = COL_ABUNDANCE) +
  geom_point(size = 2, colour = COL_ABUNDANCE) +
  antibiotic_marker_date() +
  facet_wrap(~ participant, scales = "free_x") +
  participant_date_scales() +
  labs(title = "Total normalized viral abundance over time",
       x = "Collection date", y = "Normalized abundance (size-factor adjusted)") +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(style_figure(p_viral_abundance))

# -----------------------------------------------------------------------
# 2.3. vOTU ASSEMBLY/CLUSTERING CONCORDANCE (CD-HIT-based)
# -----------------------------------------------------------------------
# NOT a persistence measure. confidence_tiers$n_samples only counts a
# vOTU as multi-sample if the SAME representative contig independently
# assembled and clustered across >=2 samples - a much stricter bar than
# "is this phage present in multiple samples." Treating this as
# persistence earlier produced a biologically implausible "near-total
# viral turnover" finding. This section is kept because tier/n_samples
# is still useful for characterizing assembly/annotation confidence
# (Aim B framing) - the real persistence comparison is Section 2.3B below.

p_cluster_concordance <- confidence_tiers |>
  count(n_samples) |>
  ggplot(aes(x = factor(n_samples), y = n)) +
  geom_col(fill = COL_VIRAL) +
  geom_text(aes(label = n), vjust = -0.3) +
  labs(title = "vOTU assembly/clustering concordance across samples (CD-HIT-based)",
       subtitle = "NOT a persistence measure - see Section 2.3B for read-mapping-based persistence",
       x = "Number of samples the SAME representative contig assembled+clustered in",
       y = "Number of vOTUs") +
  theme_thesis()

print(style_figure(p_cluster_concordance))

# Detail table for the vOTUs concordant across >=2 samples at assembly time
concordant_votu_table <- confidence_tiers |>
  filter(n_samples >= 2) |>
  select(votu_id, representative_contig, tier, tool_support, n_samples, samples, note)

print(concordant_votu_table)

# -----------------------------------------------------------------------
# 2.3B. VIRAL vs BACTERIAL PERSISTENCE - READ-MAPPING-BASED (the real
#       Aim C persistence comparison)
# -----------------------------------------------------------------------
# Uses detection_freq_combined from Part 1 (CoverM abundance > 0 for viral,
# Bracken new_est_reads > 0 for bacterial; each vOTU/taxon's detection
# count converted to a 0-1 FREQUENCY of that participant's own samples,
# since D3/G2/L1 have different total sample counts).
#
# Validated finding (re-confirmed 2026-09-25 on bacterial_long AFTER the
# phage-contamination exclusion in Part 1 Section 1.2a - numbers barely
# moved from the pre-cleanup check, which is itself good evidence the
# finding is robust to that contamination): the bacteriome carries a
# substantial stable "core" fraction the virome largely lacks. Detected
# in >=80% of a participant's samples: 40.8% of bacterial taxa vs only
# 4.3% of vOTUs - roughly a 9-10x difference. Transient rates (<=20% of
# samples) are far more similar (27.1% bacterial vs 26.3% viral), so the
# gap is specifically in the high-persistence "core" tail, not an
# inflated transient population on the viral side. Confirmed to
# replicate WITHIN each of D3, G2, and L1 individually (not a pooling
# artifact) via per-participant Kolmogorov-Smirnov tests, all p <
# 0.0001 (R notes "p-value will be approximate in the presence of ties" -
# expected with discrete detection-frequency ratios, doesn't affect the
# conclusion).

# ECDF - best for comparing overall distribution shape/shift
p_detection_freq_ecdf <- ggplot(detection_freq_combined, aes(x = detection_freq, colour = domain)) +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  facet_wrap(~ participant) +
  scale_colour_manual(values = c("Viral" = COL_VIRAL, "Bacterial" = COL_BACTERIAL)) +
  labs(title = "Detection frequency (read-mapping-based): viral vs bacterial",
       subtitle = "Proportion of that participant's own samples each vOTU/taxon was detected in",
       x = "Detection frequency", y = "Cumulative proportion", colour = NULL) +
  theme_thesis()

print(style_figure(p_detection_freq_ecdf))

# Density histogram - shows the bimodal "core + transient" bacterial
# shape (and its absence in viral) more intuitively than the ECDF
p_detection_freq_hist <- ggplot(detection_freq_combined, aes(x = detection_freq, fill = domain)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.55, bins = 10) +
  facet_wrap(~ participant) +
  scale_fill_manual(values = c("Viral" = COL_VIRAL, "Bacterial" = COL_BACTERIAL)) +
  labs(title = "Detection frequency distribution (density-normalized): viral vs bacterial",
       x = "Detection frequency", y = "Density", fill = NULL) +
  theme_thesis()

print(style_figure(p_detection_freq_hist))

# Per-participant Kolmogorov-Smirnov test - confirms viral vs bacterial
# detection-frequency distributions differ significantly WITHIN each
# individual, not just when pooled
ks_by_participant <- detection_freq_combined |>
  group_by(participant) |>
  summarise(
    ks_p = ks.test(
      detection_freq[domain == "Viral"],
      detection_freq[domain == "Bacterial"]
    )$p.value,
    .groups = "drop"
  )
print(ks_by_participant)

# Core / intermediate / transient split - the headline effect-size numbers
detection_freq_combined |>
  mutate(category = case_when(
    detection_freq <= 0.2 ~ "transient (<=20% of samples)",
    detection_freq >= 0.8 ~ "core (>=80% of samples)",
    TRUE ~ "intermediate"
  )) |>
  count(domain, category) |>
  group_by(domain) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  arrange(domain, category) |>
  print(n = Inf)

# -----------------------------------------------------------------------
# 2.4. BACTERIAL DIVERSITY OVER TIME (mirrors viral diversity plot)
# -----------------------------------------------------------------------

bacterial_diversity <- bacterial_long |>
  group_by(sample_id, participant, timepoint) |>
  summarise(
    richness = sum(new_est_reads > 0),
    shannon  = {
      p <- new_est_reads[new_est_reads > 0] / sum(new_est_reads)
      -sum(p * log(p))
    },
    .groups = "drop"
  ) |>
  left_join(sample_dates, by = "sample_id")

p_bacterial_shannon <- ggplot(bacterial_diversity, aes(x = collection_date, y = shannon)) +
  geom_line(colour = COL_BACTERIAL) +
  geom_point(size = 2, colour = COL_BACTERIAL) +
  antibiotic_marker_date() +
  facet_wrap(~ participant, scales = "free_x") +
  participant_date_scales() +
  labs(title = "Bacterial Shannon diversity over time",
       x = "Collection date", y = "Shannon diversity") +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(style_figure(p_bacterial_shannon))

# -----------------------------------------------------------------------
# 2.5. VIRAL vs BACTERIAL DIVERSITY - COMBINED VIEW
# -----------------------------------------------------------------------
# Useful for eyeballing whether viral and bacterial diversity move
# together over time, or diverge - especially interesting around L1's
# antibiotic window.

combined_diversity <- bind_rows(
  viral_diversity    |> mutate(domain = "Viral"),
  bacterial_diversity |> mutate(domain = "Bacterial")
)

p_combined <- ggplot(combined_diversity, aes(x = collection_date, y = shannon, colour = domain)) +
  geom_line() +
  geom_point(size = 2) +
  antibiotic_marker_date() +
  facet_wrap(~ participant, scales = "free_x") +
  scale_colour_manual(values = c("Viral" = COL_VIRAL, "Bacterial" = COL_BACTERIAL)) +
  participant_date_scales() +
  labs(title = "Viral vs bacterial Shannon diversity over time",
       x = "Collection date", y = "Shannon diversity", colour = NULL) +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(style_figure(p_combined))

# -----------------------------------------------------------------------
# 2.6. AIM A FORMAL TREND TEST: is the diversity/richness trend over time
#      statistically real, or just visual?
# -----------------------------------------------------------------------
# Per-participant, per-domain, per-metric Kendall trend test (cor.test()
# against timepoint, method = "kendall") - the standard dependency-free
# way to run a Mann-Kendall-style monotonic trend test on a single time
# series in base R (Kendall's tau vs time), without adding a package for
# an otherwise simple test.
#
# Kept PER-PARTICIPANT (never pooled) for the same reason as everywhere
# else in this file: the study's focus is within-individual temporal
# change, and pooling 3 people's baselines together would conflate
# between-person differences with an actual time trend. This also avoids
# needing a mixed model with only 3 participant levels, which is too few
# for a random-intercept variance component to be estimated sensibly.
#
# Unlike Part 5's co-occurrence heatmaps (explicitly exploratory,
# uncorrected), this IS meant as a formal claim ("is the trend real"), so
# a Benjamini-Hochberg correction is applied across all 12 tests run here
# (3 participants x 2 domains x 2 metrics) to control the false discovery
# rate.

diversity_trend_tests <- combined_diversity |>
  pivot_longer(cols = c(richness, shannon), names_to = "metric", values_to = "value") |>
  group_by(participant, domain, metric) |>
  summarise(
    mk = list(suppressWarnings(cor.test(timepoint, value, method = "kendall"))),
    n_timepoints = n(),
    .groups = "drop"
  ) |>
  mutate(
    kendall_tau = purrr::map_dbl(mk, ~ unname(.x$estimate)),
    kendall_p   = purrr::map_dbl(mk, ~ .x$p.value)
  ) |>
  select(-mk) |>
  mutate(kendall_p_adj_BH = p.adjust(kendall_p, method = "BH")) |>
  arrange(domain, metric, participant)

message("=== Aim A: per-participant Kendall trend test (richness/Shannon vs timepoint), BH-corrected across all 12 tests ===")
print(diversity_trend_tests, n = Inf)
message("NOTE: n_timepoints per test is only 8-10 (D3), 8 (G2), or 9 (L1) - Kendall's tau is a nonparametric rank test appropriate for this small n and doesn't assume linearity, but a 'significant' trend here should still be read as a monotonic tendency over a short series, not a precise rate of change. With ties present (repeated exact values), cor.test() falls back to a normal approximation rather than an exact p-value - same caveat already noted for ks.test() above.")


# #############################################################################
# PART 3: TAXONOMY
# (originally taxonomy.R)
# Viral (PhaBOX2 Lineage) and bacterial (Bracken, genus parsed from species
# names) taxonomic composition. Depends on Part 1's viral_long, taxonomy,
# bacterial_long, AND Part 2's antibiotic_marker() helper (used in Sections
# 3.5B, 3.6, 3.7 below).
# #############################################################################

# -----------------------------------------------------------------------
# 3.1. VIRAL - classification DEPTH (non-temporal: taxonomy is a fixed
#      per-vOTU property, like lifestyle in aim_b_functional.R)
# -----------------------------------------------------------------------
# PhaBOX2's Genus column reads "-" for most vOTUs, but the
# underlying Lineage string often still carries a real classification at
# a coarser rank (Family, Order, Class, Phylum...) even when Genus is
# unresolved - e.g. some vOTUs are classified to Family (Autotranscript-
# aviridae) or Class (Caudoviricetes) despite showing Genus = "-". Using
# Genus alone would misreport those as "unclassified" when they aren't.
# This instead extracts the DEEPEST rank each vOTU was actually assigned
# to, in the order PhaBOX2's Lineage strings use:
#   root/superkingdom < realm < kingdom < phylum < class < order <
#   family < subfamily < genus < species
# (Lineage strings use "superkingdom" or "acellular root" inconsistently
# for the same top rank - both are the top of the hierarchy for a virus
# lineage, so both get normalized to "root" here.)

rank_order <- c("unclassified", "root", "realm", "kingdom", "phylum",
                "class", "order", "family", "subfamily", "genus", "species")

extract_deepest_rank <- function(lineage) {
  if (is.na(lineage) || lineage == "-") return("unclassified")
  parts <- str_split(lineage, ";")[[1]]
  last_rank <- str_extract(parts[length(parts)], "^[^:]+")
  if (last_rank %in% c("acellular root", "superkingdom")) return("root")
  last_rank
}

# Companion to extract_deepest_rank: the ACTUAL NAME at the deepest rank
# reached (e.g. "Family: Autotranscriptaviridae"), not just which rank
# level was reached. Needed for anything that labels an individual vOTU
# (e.g. Section 3.7's line-chart panels) - "reached family rank" on its
# own says nothing about which virus this actually is.
extract_deepest_taxon_label <- function(lineage) {
  if (is.na(lineage) || lineage == "-") return("Unclassified")
  parts <- str_split(lineage, ";")[[1]]
  last_part <- parts[length(parts)]
  rank_name <- str_extract(last_part, "^[^:]+")
  value     <- str_extract(last_part, "(?<=:).+")
  if (rank_name %in% c("acellular root", "superkingdom")) {
    return("Unclassified")  # root-only ("Viruses") carries no real information
  }
  paste0(str_to_title(rank_name), ": ", value)
}

taxonomy <- taxonomy |>
  mutate(
    deepest_rank  = map_chr(Lineage, extract_deepest_rank),
    deepest_rank  = factor(deepest_rank, levels = rank_order),
    taxon_label   = map_chr(Lineage, extract_deepest_taxon_label)
  )

message("=== Viral vOTU classification depth (deepest rank reached) ===")
taxonomy |>
  count(deepest_rank) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  print(n = Inf)

p_viral_classification_depth <- taxonomy |>
  count(deepest_rank) |>
  ggplot(aes(x = deepest_rank, y = n)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = n), vjust = -0.3) +
  labs(title = "vOTU classification depth (PhaBOX2 / PhaGCN)",
       subtitle = sprintf("Deepest taxonomic rank each of the %d vOTUs could be confidently assigned to", nrow(taxonomy)),
       x = "Deepest rank reached", y = "Number of vOTUs") +
  theme_thesis()

print(style_figure(p_viral_classification_depth))

# Join deepest_rank and taxon_label back into viral_long for anything
# downstream that wants them (e.g. host-association work later in Aim C,
# or Section 3.7's line-chart panel labels)
viral_long <- viral_long |>
  left_join(taxonomy |> select(contig, deepest_rank, taxon_label), by = "contig")

# -----------------------------------------------------------------------
# 3.1B. VIRAL - classification depth PER PARTICIPANT
# -----------------------------------------------------------------------
# IMPORTANT: "per participant" here means which vOTUs that participant
# actually HAS, determined by read-mapping detection (abundance > 0 in
# >=1 of their own samples) - the same logic established for persistence
# (Part 1 Section 1.4). It is NOT which participant's sample happened to
# assemble the representative contig - a vOTU can be detected via
# mapping in a participant's samples even when its representative
# assembled from someone else's sample. A vOTU detected in more than one
# participant will appear in more than one participant's breakdown
# below, which is correct, not a double-count error.

viral_detected_by_participant <- viral_long |>
  filter(abundance > 0) |>
  distinct(contig, participant, deepest_rank, Genus)

message("=== Viral vOTU classification depth, per participant (mapping-based) ===")
viral_detected_by_participant |>
  count(participant, deepest_rank) |>
  group_by(participant) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  arrange(participant, deepest_rank) |>
  print(n = Inf)

p_viral_depth_by_participant <- viral_detected_by_participant |>
  count(participant, deepest_rank) |>
  ggplot(aes(x = deepest_rank, y = n, fill = participant)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = PARTICIPANT_COLORS) +
  labs(title = "vOTU classification depth per participant (mapping-based detection)",
       x = "Deepest rank reached", y = "Number of vOTUs detected", fill = NULL) +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(style_figure(p_viral_depth_by_participant))

# -----------------------------------------------------------------------
# 3.2. VIRAL - genus composition, classified vOTUs only (non-temporal)
# -----------------------------------------------------------------------
# Only a small minority of vOTUs reach genus level - a table is more honest
# here than a bar chart, which would imply meaningful category sizes at
# n=1-4 per genus.

n_genus_called <- sum(taxonomy$Genus != "-")
message(sprintf("=== Classified vOTU genera (n = %d/%d, %.1f%%) ===",
                n_genus_called, nrow(taxonomy), 100 * n_genus_called / nrow(taxonomy)))
taxonomy |>
  filter(Genus != "-") |>
  count(Genus, sort = TRUE) |>
  print(n = Inf)

message("=== Classified vOTU genera detected per participant (mapping-based) ===")
viral_detected_by_participant |>
  filter(Genus != "-") |>
  count(participant, Genus, sort = TRUE) |>
  print(n = Inf)

# -----------------------------------------------------------------------
# 3.3. BACTERIAL - genus parsing
# -----------------------------------------------------------------------
# Bracken only reports at species level (573 distinct taxa after the
# phage-contamination exclusion in Part 1 Section 1.2a, all
# taxonomy_lvl == "S") - no Phylum/Genus column exists directly, so
# Genus is parsed as the first word of the species name (standard
# binomial convention). Checked for parsing risk first:
#   - "[Ruminococcus] torques" etc. (6 taxa): square brackets denote a
#     genus in taxonomic flux, not a different genus - stripped so
#     these merge correctly with any unbracketed entries of the same
#     genus.
#   - "Candidatus Nanosynbacter..." etc. (4 taxa): a naive first-word
#     split would bucket these under the literal word "Candidatus"
#     instead of the real genus - prefix stripped first.
# (Phage/viral contamination - e.g. "uncultured phage WW-nAnB", and 27
# properly-named ICTV phage genera like "Endlipuvirus intestinihominis" -
# is now excluded upstream in Part 1 Section 1.2a, not here, so it's
# cleaned before ANY downstream use of bacterial_long, including the
# persistence comparison in Part 1 Section 1.4.)

bacterial_long <- bacterial_long |>
  mutate(
    bacterial_genus = bacterial_taxon |>
      str_remove("^Candidatus ") |>
      str_remove_all("[\\[\\]]") |>
      word(1)
  )

message(sprintf(
  "Bacterial genera: %d distinct (parsed from %d distinct species-level taxa)",
  n_distinct(bacterial_long$bacterial_genus),
  n_distinct(bacterial_long$bacterial_taxon)
))

# -----------------------------------------------------------------------
# 3.4. BACTERIAL - overall genus composition (non-temporal, pooled across
#      all samples) - top 15 genera by total abundance + "Other"
# -----------------------------------------------------------------------

top_genera <- bacterial_long |>
  group_by(bacterial_genus) |>
  summarise(total_reads = sum(new_est_reads), .groups = "drop") |>
  slice_max(total_reads, n = 15) |>
  pull(bacterial_genus)

bacterial_long <- bacterial_long |>
  mutate(genus_grouped = if_else(bacterial_genus %in% top_genera, bacterial_genus, "Other"))

p_bacterial_genus_overall <- bacterial_long |>
  group_by(genus_grouped) |>
  summarise(total_reads = sum(new_est_reads), .groups = "drop") |>
  mutate(genus_grouped = fct_reorder(genus_grouped, total_reads)) |>
  ggplot(aes(x = genus_grouped, y = total_reads)) +
  geom_col(fill = COL_BACTERIAL) +
  coord_flip() +
  labs(title = "Overall bacterial genus composition (pooled across all samples)",
       subtitle = "Top 15 genera by total Bracken read count; remainder grouped as 'Other'",
       x = NULL, y = "Total estimated reads") +
  theme_thesis()

print(style_figure(p_bacterial_genus_overall))

# -----------------------------------------------------------------------
# 3.4B. BACTERIAL - overall genus composition PER PARTICIPANT (non-temporal,
#       pooled across each participant's own timepoints)
# -----------------------------------------------------------------------
# Uses the SAME top-15-genera-overall grouping as Section 3.4 (not a
# per-participant top-15) so the three bars share identical categories -
# necessary for a fair visual comparison across participants.

p_bacterial_genus_by_participant <- bacterial_long |>
  group_by(participant, genus_grouped) |>
  summarise(total_reads = sum(new_est_reads), .groups = "drop") |>
  ggplot(aes(x = participant, y = total_reads, fill = genus_grouped)) +
  geom_col(position = "fill") +
  scale_fill_viridis_d(option = "viridis") +
  labs(title = "Bacterial genus composition per participant",
       subtitle = "Pooled across each participant's own timepoints; top 15 genera overall + 'Other'",
       x = NULL, y = "Proportion of reads", fill = NULL) +
  theme_thesis() +
  theme(legend.position = "bottom")

print(style_figure(p_bacterial_genus_by_participant))

# -----------------------------------------------------------------------
# 3.5. BACTERIAL - genus composition OVER TIME (abundance-weighted,
#      faceted by participant) - mirrors the viral functional-composition
#      figure in aim_b_functional.R Section 3. Already per-participant via
#      the facet, unlike Section 3.4B above which pools each participant's
#      timepoints into a single bar.
# -----------------------------------------------------------------------

genus_over_time <- bacterial_long |>
  group_by(participant, sample_id, timepoint, genus_grouped) |>
  summarise(summed_reads = sum(new_est_reads), .groups = "drop")

# NOTE ON X-AXIS: this and Section 3.5B stay on `timepoint` (sample order)
# rather than moving to collection_date like the trend/line plots elsewhere
# in this file. These are stacked-proportion snapshots (one full-height bar
# per sample) - on a real calendar axis, L1's uneven 14-49 day gaps between
# visits would force very unequal bar widths/gaps, which reads as "some
# samples matter more" rather than reflecting anything about the
# composition itself. Sample order is the more honest axis for this
# specific chart type; integer_breaks() (Part 2) keeps the tick labels as
# whole sample numbers rather than ggplot's default fractional breaks.
p_bacterial_genus_over_time <- ggplot(genus_over_time,
                                      aes(x = timepoint, y = summed_reads, fill = genus_grouped)) +
  geom_col(position = "fill") +
  antibiotic_marker() +
  facet_wrap(~ participant, scales = "free_x") +
  scale_x_continuous(breaks = integer_breaks) +
  scale_fill_viridis_d(option = "viridis") +
  labs(title = "Bacterial genus composition over time (top 15 genera + Other)",
       x = "Timepoint (sample order)", y = "Proportion of reads", fill = NULL) +
  theme_thesis() +
  theme(legend.position = "bottom")

print(style_figure(p_bacterial_genus_over_time))

# -----------------------------------------------------------------------
# 3.5B. BACTERIAL - genus composition over time, ONE FIGURE PER PARTICIPANT
# -----------------------------------------------------------------------
# Same data/categories as Section 3.5, but as three separate full-size
# plots instead of small facets - easier to read with 16 stacked
# categories per panel. IMPORTANT: antibiotic_marker() is built to work
# WITHIN a facet (ggplot only draws it in the L1 panel because its data
# carries a participant column matched against the facet). Outside a
# facet, on a single-participant plot, that matching doesn't happen -
# it would draw on all three plots regardless of participant - so it's
# applied conditionally here, only when plotting L1.

for (p in unique(genus_over_time$participant)) {
  
  plot_data <- genus_over_time |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = timepoint, y = summed_reads, fill = genus_grouped)) +
    geom_col(position = "fill") +
    scale_x_continuous(breaks = integer_breaks) +
    scale_fill_viridis_d(option = "viridis") +
    labs(title = paste0("Bacterial genus composition over time - ", p),
         subtitle = "Top 15 genera overall + 'Other' (proportion of reads, abundance-weighted)",
         x = "Timepoint (sample order)", y = "Proportion of reads", fill = NULL) +
    theme_thesis() +
    theme(legend.position = "bottom")
  
  if (p == "L1") gg <- gg + antibiotic_marker()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 3.6. BACTERIAL - trajectories of the most variable genera, per participant
# -----------------------------------------------------------------------
# "Biologically interesting / unexpectedly changing" is operationalized
# as: highest coefficient of variation (CV = sd/mean) among genera
# detected in >=4 of that participant's own timepoints (the floor avoids
# a genus that's a one-off spike/near-zero-everywhere dominating the
# "most variable" list purely by being noisy at low abundance).
#
# CV=0 FIX (2026-09-29): bacterial_long only carries a row where Bracken
# actually detected a genus in a sample - there's no row for samples
# where it wasn't. Computing mean/sd straight off bacterial_long (the
# original approach) therefore averaged only over the DETECTED samples,
# which produced a mechanical CV=0 floor whenever a genus was detected
# at the exact same low count every time it appeared (e.g. always
# exactly 1 read) - that reads as "perfectly stable" but is really just
# an artefact of ignoring every sample where the genus was actually
# absent. Fixed below by explicitly zero-filling: every genus a
# participant ever had gets one row per sample of THAT participant
# (reads = 0 where undetected), so mean/sd/cv are computed across the
# participant's full sample set, not just the detected subset.

genus_reads_detected <- bacterial_long |>
  group_by(participant, bacterial_genus, sample_id, timepoint) |>
  summarise(reads = sum(new_est_reads), .groups = "drop")

participant_samples_all <- bacterial_long |>
  distinct(participant, sample_id, timepoint)

genus_reads_full <- genus_reads_detected |>
  distinct(participant, bacterial_genus) |>
  left_join(participant_samples_all, by = "participant", relationship = "many-to-many") |>
  left_join(genus_reads_detected, by = c("participant", "bacterial_genus", "sample_id", "timepoint")) |>
  mutate(reads = replace_na(reads, 0)) |>
  left_join(sample_dates, by = "sample_id")

genus_variability <- genus_reads_full |>
  group_by(participant, bacterial_genus) |>
  summarise(
    total_samples  = n(),
    n_detected     = sum(reads > 0),
    detection_freq = round(n_detected / total_samples, 2),
    mean_reads     = mean(reads),
    sd_reads       = sd(reads),
    cv             = sd_reads / mean_reads,
    .groups = "drop"
  )

candidate_genera <- genus_variability |>
  group_by(participant) |>
  filter(n_detected >= 4) |>
  slice_max(cv, n = 8) |>
  arrange(participant, desc(cv))

message("=== Most variable bacterial genera per participant (zero-filled CV; candidates for the line charts below) ===")
print(candidate_genera, n = Inf)

# Each genus gets its OWN y-axis (facet_wrap scales="free_y") rather
# than sharing one axis or a log transform - genus abundances here span
# several orders of magnitude within a single person (e.g. D3's
# Bifidobacterium mean ~3645 vs Lactiplantibacillus mean ~4.5), so a
# shared axis would flatten the low-abundance genera to invisible flat
# lines. This does mean the panels aren't magnitude-comparable to each
# other - only fair, since these genera's typical abundance level isn't
# a meaningful thing to compare directly anyway; what matters here is
# each genus's own trend shape over time. Now plotted from the
# zero-filled genus_reads_full, so a genus's true absence at a given
# timepoint shows as a real dip to 0 rather than a silently skipped gap.

genus_line_data <- genus_reads_full |>
  semi_join(candidate_genera, by = c("participant", "bacterial_genus"))

for (p in unique(genus_line_data$participant)) {
  
  plot_data <- genus_line_data |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = reads)) +
    geom_line(colour = COL_BACTERIAL_DARK) +
    geom_point(size = 1.8, colour = COL_BACTERIAL_DARK) +
    facet_wrap(~ bacterial_genus, scales = "free_y") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Most variable bacterial genera over time - ", p),
         subtitle = "Top 8 genera by coefficient of variation (zero-filled; detected in >=4 of this participant's own timepoints)",
         x = "Collection date", y = "Estimated reads (Bracken)") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold.italic"))
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 3.6B. BACTERIAL - trajectories of the LEAST variable (most stable)
#       genera, per participant - the mirror image of Section 3.6.
# -----------------------------------------------------------------------
# Same zero-filled CV (sd/mean) and >=4-detection floor as Section 3.6
# (reuses genus_variability directly), just slice_min() instead of
# slice_max(), so the two lists are directly comparable. detection_freq
# is still reported alongside CV: with CV now computed on the full
# sample set, a low CV and a low detection_freq can no longer coexist
# purely as an artefact the way they could before (a genus mostly absent
# now pulls its own mean down and its CV up, per the fix above), but
# detection_freq remains useful context for how much data actually
# supports the "stable" call for that genus.

candidate_genera_stable <- genus_variability |>
  group_by(participant) |>
  filter(n_detected >= 4) |>
  slice_min(cv, n = 8) |>
  arrange(participant, cv)

message("=== Least variable (most stable) bacterial genera per participant (zero-filled CV; candidates for the line charts below) ===")
print(candidate_genera_stable, n = Inf)

genus_line_data_stable <- genus_reads_full |>
  semi_join(candidate_genera_stable, by = c("participant", "bacterial_genus"))

for (p in unique(genus_line_data_stable$participant)) {
  
  plot_data <- genus_line_data_stable |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = reads)) +
    geom_line(colour = COL_BACTERIAL_LIGHT) +
    geom_point(size = 1.8, colour = COL_BACTERIAL_LIGHT) +
    facet_wrap(~ bacterial_genus, scales = "free_y") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Least variable (most stable) bacterial genera over time - ", p),
         subtitle = "Bottom 8 genera by coefficient of variation (zero-filled; detected in >=4 of this participant's own timepoints)",
         x = "Collection date", y = "Estimated reads (Bracken)") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold.italic"))
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 3.6C. BACTERIAL - trajectories of the most variable SPECIES, per
#       participant - species-level counterpart to Section 3.6. Same
#       zero-filled CV logic, grouped on bacterial_taxon (the raw
#       Bracken species-level name) instead of the genus parsed from it.
# -----------------------------------------------------------------------

species_reads_detected <- bacterial_long |>
  group_by(participant, bacterial_taxon, sample_id, timepoint) |>
  summarise(reads = sum(new_est_reads), .groups = "drop")

species_reads_full <- species_reads_detected |>
  distinct(participant, bacterial_taxon) |>
  left_join(participant_samples_all, by = "participant", relationship = "many-to-many") |>
  left_join(species_reads_detected, by = c("participant", "bacterial_taxon", "sample_id", "timepoint")) |>
  mutate(reads = replace_na(reads, 0)) |>
  left_join(sample_dates, by = "sample_id")

species_variability <- species_reads_full |>
  group_by(participant, bacterial_taxon) |>
  summarise(
    total_samples  = n(),
    n_detected     = sum(reads > 0),
    detection_freq = round(n_detected / total_samples, 2),
    mean_reads     = mean(reads),
    sd_reads       = sd(reads),
    cv             = sd_reads / mean_reads,
    .groups = "drop"
  )

candidate_species <- species_variability |>
  group_by(participant) |>
  filter(n_detected >= 4) |>
  slice_max(cv, n = 8) |>
  arrange(participant, desc(cv))

message("=== Most variable bacterial SPECIES per participant (zero-filled CV; candidates for the line charts below) ===")
print(candidate_species, n = Inf)

species_line_data <- species_reads_full |>
  semi_join(candidate_species, by = c("participant", "bacterial_taxon"))

for (p in unique(species_line_data$participant)) {
  
  plot_data <- species_line_data |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = reads)) +
    geom_line(colour = COL_BACTERIAL_DARK) +
    geom_point(size = 1.8, colour = COL_BACTERIAL_DARK) +
    facet_wrap(~ bacterial_taxon, scales = "free_y") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Most variable bacterial species over time - ", p),
         subtitle = "Top 8 species by coefficient of variation (zero-filled; detected in >=4 of this participant's own timepoints)",
         x = "Collection date", y = "Estimated reads (Bracken)") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 7, face = "bold.italic"))
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 3.6D. BACTERIAL - trajectories of the LEAST variable (most stable)
#       SPECIES, per participant - species-level counterpart to 3.6B.
# -----------------------------------------------------------------------

candidate_species_stable <- species_variability |>
  group_by(participant) |>
  filter(n_detected >= 4) |>
  slice_min(cv, n = 8) |>
  arrange(participant, cv)

message("=== Least variable (most stable) bacterial SPECIES per participant (zero-filled CV; candidates for the line charts below) ===")
print(candidate_species_stable, n = Inf)

species_line_data_stable <- species_reads_full |>
  semi_join(candidate_species_stable, by = c("participant", "bacterial_taxon"))

for (p in unique(species_line_data_stable$participant)) {
  
  plot_data <- species_line_data_stable |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = reads)) +
    geom_line(colour = COL_BACTERIAL_LIGHT) +
    geom_point(size = 1.8, colour = COL_BACTERIAL_LIGHT) +
    facet_wrap(~ bacterial_taxon, scales = "free_y") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Least variable (most stable) bacterial species over time - ", p),
         subtitle = "Bottom 8 species by coefficient of variation (zero-filled; detected in >=4 of this participant's own timepoints)",
         x = "Collection date", y = "Estimated reads (Bracken)") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 7, face = "bold.italic"))
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 3.7. VIRAL - trajectories of the most variable vOTUs, per participant
# -----------------------------------------------------------------------
# Same CV-based selection logic as Section 3.6, applied to read-mapping
# abundance (viral_long$abundance, CoverM trimmed-mean coverage) instead
# of Bracken reads. A vOTU appearing in more than one participant's
# figure here (e.g. a vOTU detected via mapping in both D3 and G2) is
# real cross-participant detection, not a duplicate/error - see the
# note in Section 3.1B.
#
# ZERO-FILL CHECK (2026-09-29, same audit that fixed the bacterial CV=0
# artifact above): unlike bacterial_long, viral_long already carries a
# zero-filled row for every contig in every sample - it's built by
# pivoting the WIDE CoverM abundance matrix long (Part 1, Section 1.1a),
# which has an explicit 0 entry for every contig not detected in a given
# sample, rather than only listing detections. Part 7 Section 7.1's bug
# fix independently confirms the same fact. So mean_abund/sd_abund below
# were already being computed across each participant's FULL sample set,
# not just the detected subset - the mechanical CV=0 floor that affected
# the bacterial genus table does not apply here. Kept as an explicit,
# separately-named computation (not silently reused from before) so
# that's verified rather than assumed.

votu_reads_full <- viral_long |>
  group_by(participant, contig, sample_id, timepoint, collection_date) |>
  summarise(abund = sum(abundance), .groups = "drop")

votu_variability <- votu_reads_full |>
  group_by(participant, contig) |>
  summarise(
    total_samples  = n(),
    n_detected     = sum(abund > 0),
    detection_freq = round(n_detected / total_samples, 2),
    mean_abund     = mean(abund),
    sd_abund       = sd(abund),
    cv             = sd_abund / mean_abund,
    .groups = "drop"
  )

candidate_votus <- votu_variability |>
  group_by(participant) |>
  filter(n_detected >= 4) |>
  slice_max(cv, n = 8) |>
  arrange(participant, desc(cv))

message("=== Most variable vOTUs per participant (candidates for the line charts below) ===")
print(candidate_votus, n = Inf)

# Lookup table: what do we actually know taxonomically about these 24
# candidates? Most will likely read "Unclassified" - PhaBOX2 only
# resolved SOME rank for a minority of the 151-vOTU catalogue overall
# (Section 3.1) - but this makes clear exactly which ones and to what
# depth, rather than silently plotting bare contig IDs.
message("=== Taxonomic classification of the selected 'most variable' vOTUs ===")
candidate_votus |>
  left_join(taxonomy |> select(contig, taxon_label), by = "contig") |>
  select(participant, contig, taxon_label, cv, n_detected) |>
  arrange(participant, desc(cv)) |>
  print(n = Inf)

votu_line_data <- votu_reads_full |>
  semi_join(candidate_votus, by = c("participant", "contig")) |>
  left_join(taxonomy |> select(contig, taxon_label), by = "contig") |>
  left_join(votu_display_lookup |> select(contig, votu_display), by = "contig") |>
  mutate(facet_label = paste0(votu_display, "\n", taxon_label))

for (p in unique(votu_line_data$participant)) {
  
  plot_data <- votu_line_data |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = abund)) +
    geom_line(colour = COL_VIRAL_DARK) +
    geom_point(size = 1.8, colour = COL_VIRAL_DARK) +
    facet_wrap(~ facet_label, scales = "free_y") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Most variable vOTUs over time - ", p),
         subtitle = "Top 8 vOTUs by coefficient of variation (detected in >=4 of this participant's own timepoints, via read-mapping)",
         x = "Collection date", y = "CoverM trimmed-mean coverage") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 3.7B. VIRAL - trajectories of the LEAST variable (most stable) vOTUs,
#       per participant - the mirror image of Section 3.7.
# -----------------------------------------------------------------------

candidate_votus_stable <- votu_variability |>
  group_by(participant) |>
  filter(n_detected >= 4) |>
  slice_min(cv, n = 8) |>
  arrange(participant, cv)

message("=== Least variable (most stable) vOTUs per participant (candidates for the line charts below) ===")
print(candidate_votus_stable, n = Inf)

message("=== Taxonomic classification of the selected 'least variable' vOTUs ===")
candidate_votus_stable |>
  left_join(taxonomy |> select(contig, taxon_label), by = "contig") |>
  select(participant, contig, taxon_label, cv, n_detected, detection_freq) |>
  arrange(participant, cv) |>
  print(n = Inf)

votu_line_data_stable <- votu_reads_full |>
  semi_join(candidate_votus_stable, by = c("participant", "contig")) |>
  left_join(taxonomy |> select(contig, taxon_label), by = "contig") |>
  left_join(votu_display_lookup |> select(contig, votu_display), by = "contig") |>
  mutate(facet_label = paste0(votu_display, "\n", taxon_label))

for (p in unique(votu_line_data_stable$participant)) {
  
  plot_data <- votu_line_data_stable |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = abund)) +
    geom_line(colour = COL_VIRAL_LIGHT) +
    geom_point(size = 1.8, colour = COL_VIRAL_LIGHT) +
    facet_wrap(~ facet_label, scales = "free_y") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Least variable (most stable) vOTUs over time - ", p),
         subtitle = "Bottom 8 vOTUs by coefficient of variation (detected in >=4 of this participant's own timepoints, via read-mapping)",
         x = "Collection date", y = "CoverM trimmed-mean coverage") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}


# #############################################################################
# PART 4: FUNCTIONAL ANNOTATION
# (originally aim_b_functional.R)
# Functional annotation visualisation (Aim B), built from Phold's per-gene
# predictions. Depends on Part 1's viral_long/confidence_tiers/BASE, and
# Part 2's antibiotic_marker().
# #############################################################################

# -----------------------------------------------------------------------
# 4.0. Load Phold's per-gene annotations
# -----------------------------------------------------------------------
phold_annotations <- read_tsv(
  file.path(BASE, "results/annotation/phold/phold_annotations.tsv"),
  show_col_types = FALSE
)

# -----------------------------------------------------------------------
# IMPORTANT STRUCTURAL NOTE:
# phold_annotations$contig_id is a vOTU REPRESENTATIVE contig, tied to
# whichever single sample that representative sequence happened to
# assemble in. Gene content is a fixed property of that one genome - it
# does NOT vary by sample. What varies by sample/timepoint is the vOTU's
# ABUNDANCE (already in viral_long). So any "functional composition over
# time" plot must weight function-category membership by real per-sample
# vOTU abundance, not facet by the representative contig's origin
# sample - the latter would misleadingly suggest a function only
# "existed" wherever its representative happened to assemble.
# -----------------------------------------------------------------------

# Map each gene's representative contig to its vOTU id
phold_votu <- phold_annotations |>
  left_join(
    confidence_tiers |> select(votu_id, representative_contig),
    by = c("contig_id" = "representative_contig")
  )

# -----------------------------------------------------------------------
# 4.1. OVERALL FUNCTIONAL COMPOSITION (non-temporal - genuinely fine as-is,
#      since this just describes gene content pooled across all vOTUs,
#      with no time claim attached)
# -----------------------------------------------------------------------

p_function_overall <- ggplot(phold_annotations, aes(x = fct_infreq(`function`))) +
  geom_bar(fill = COL_VIRAL) +
  coord_flip() +
  labs(title = "Overall functional category composition (all genes, all vOTUs)",
       x = NULL, y = "Number of genes") +
  theme_thesis()

print(style_figure(p_function_overall))

# -----------------------------------------------------------------------
# 4.2. FUNCTIONAL COMPOSITION BY LIFESTYLE (non-temporal - lifestyle is
#      also a fixed per-vOTU property, no time-weighting needed)
# -----------------------------------------------------------------------

votu_lifestyle <- viral_long |>
  distinct(contig, TYPE)

phold_votu_lifestyle <- phold_votu |>
  left_join(votu_lifestyle, by = c("contig_id" = "contig")) |>
  filter(TYPE != "-")

p_function_by_lifestyle <- ggplot(phold_votu_lifestyle,
                                  aes(x = `function`, fill = TYPE)) +
  geom_bar(position = "fill") +
  coord_flip() +
  scale_fill_viridis_d(option = "viridis") +
  labs(title = "Functional category composition by phage lifestyle",
       x = NULL, y = "Proportion", fill = "Lifestyle") +
  theme_thesis()

print(style_figure(p_function_by_lifestyle))

# -----------------------------------------------------------------------
# 4.3. FUNCTIONAL COMPOSITION OVER TIME (abundance-weighted - the corrected
#      temporal version)
# -----------------------------------------------------------------------
# For each vOTU, which function categories does it carry at least one
# gene of? ("unknown function" excluded here since it carries no
# interpretable signal - remove this filter if you want it included.)

votu_function_presence <- phold_votu |>
  distinct(votu_id, `function`) |>
  filter(`function` != "unknown function")

# Expand each vOTU's real per-sample abundance across every function
# category it carries, then sum abundance per sample/timepoint/category.
# NOTE: viral_long already carries votu_id (joined in wholesale from
# confidence_tiers back in Part 1, Section 1.1d) - no need to rejoin
# confidence_tiers here. An earlier version of this block did rejoin it,
# which silently collided with the existing votu_id column (dplyr renamed
# both to votu_id.x/votu_id.y instead of erroring at that join) and only
# surfaced as a failure on the NEXT join, which expected a plain "votu_id"
# column that no longer existed.
functional_abundance_long <- viral_long |>
  inner_join(votu_function_presence, by = "votu_id", relationship = "many-to-many") |>
  filter(abundance > 0)

functional_summary <- functional_abundance_long |>
  group_by(participant, sample_id, timepoint, `function`) |>
  summarise(summed_abundance = sum(abundance), .groups = "drop") |>
  left_join(sample_dates, by = "sample_id")

# X-AXIS NOTE: stays on `timepoint` (sample order), same reasoning as the
# bacterial genus composition plots in Part 3 (Section 3.5) - this is a
# stacked-proportion snapshot per sample, and real calendar spacing would
# distort bar widths/gaps rather than clarify anything. integer_breaks()
# keeps tick labels as whole sample numbers.
p_function_over_time <- ggplot(functional_summary,
                               aes(x = timepoint, y = summed_abundance, fill = `function`)) +
  geom_col(position = "fill") +
  antibiotic_marker() +
  facet_wrap(~ participant, scales = "free_x") +
  scale_x_continuous(breaks = integer_breaks) +
  scale_fill_viridis_d(option = "viridis") +
  labs(title = "Functional category composition over time (abundance-weighted)",
       x = "Timepoint (sample order)", y = "Proportion of abundance", fill = NULL) +
  theme_thesis() +
  theme(legend.position = "bottom")

print(style_figure(p_function_over_time))

# -----------------------------------------------------------------------
# 4.4. AUXILIARY METABOLIC GENES (AMGs) SPECIFICALLY, OVER TIME
# -----------------------------------------------------------------------
# AMGs (phage-carried genes that directly affect host bacterial
# metabolism) are often of particular interest - a dedicated trend line,
# with L1's antibiotic window marked, since AMG carriage shifting around
# a perturbation would be a genuinely interesting finding.

amg_over_time <- functional_summary |>
  filter(`function` == "moron, auxiliary metabolic gene and host takeover")

p_amg_over_time <- ggplot(amg_over_time, aes(x = collection_date, y = summed_abundance)) +
  geom_line(colour = COL_AMG) +
  geom_point(size = 2, colour = COL_AMG) +
  antibiotic_marker_date() +
  facet_wrap(~ participant, scales = "free_x") +
  participant_date_scales() +
  labs(title = "Auxiliary metabolic gene (AMG)-carrying vOTU abundance over time",
       x = "Collection date", y = "Summed abundance (AMG-carrying vOTUs)") +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(style_figure(p_amg_over_time))

# -----------------------------------------------------------------------
# 4.5. AIM B FORMAL TEST: does lifestyle significantly predict functional
#      category distribution?
# -----------------------------------------------------------------------
# Unit-of-replication note: p_function_by_lifestyle (Section 4.2) plots
# GENE-level proportions, which is fine for a descriptive picture but
# would be pseudoreplicated as a hypothesis test - a single vOTU carrying
# many genes of the same function contributes many correlated "votes" to
# one cell of the contingency table, inflating the effective sample size
# and the significance of any association. This test instead uses
# vOTU-level PRESENCE (does this vOTU carry >=1 gene of this function,
# yes/no), reusing Section 4.3's votu_function_presence - one row per
# (votu, function) actually carried - joined to lifestyle at the vOTU
# level (not the gene level) so each vOTU contributes to a cell at most
# once per function category.

votu_lifestyle_by_votu <- confidence_tiers |>
  select(votu_id, representative_contig) |>
  left_join(votu_lifestyle, by = c("representative_contig" = "contig")) |>
  filter(!is.na(TYPE), TYPE != "-") |>
  distinct(votu_id, TYPE)

votu_function_lifestyle <- votu_function_presence |>
  inner_join(votu_lifestyle_by_votu, by = "votu_id")

lifestyle_function_table <- table(votu_function_lifestyle$`function`, votu_function_lifestyle$TYPE)

message(sprintf("=== Lifestyle x function contingency table (vOTU-level; %d vOTUs with both a lifestyle call and >=1 annotated function) ===",
                n_distinct(votu_function_lifestyle$votu_id)))
print(lifestyle_function_table)
# One Fisher's exact test per function category (2x2: vOTU carries the
# function yes/no x temperate/virulent). A vOTU carries several functions,
# so a single chi-square on the function x lifestyle table above would count
# the same vOTU in up to nine rows (non-independent observations,
# anti-conservative). Here each test uses every lifestyle-called vOTU
# exactly once; the nine p-values are BH-adjusted. The nine tests are not
# independent of each other either (same vOTUs), so BH is a convenience
# correction, not an exact one.
lifestyle_votu_n <- votu_lifestyle_by_votu |>
  filter(votu_id %in% votu_function_lifestyle$votu_id) |>
  count(TYPE, name = "n_total")
n_temp_votus <- lifestyle_votu_n$n_total[lifestyle_votu_n$TYPE == "temperate"]
n_vir_votus  <- lifestyle_votu_n$n_total[lifestyle_votu_n$TYPE == "virulent"]

lifestyle_function_fisher <- votu_function_lifestyle |>
  distinct(votu_id, `function`, TYPE) |>
  count(`function`, TYPE, name = "n_with") |>
  tidyr::pivot_wider(names_from = TYPE, values_from = n_with, values_fill = 0) |>
  mutate(
    temperate_with = temperate, virulent_with = virulent,
    temperate_without = n_temp_votus - temperate_with,
    virulent_without  = n_vir_votus  - virulent_with,
    prop_temperate = temperate_with / n_temp_votus,
    prop_virulent  = virulent_with  / n_vir_votus
  ) |>
  rowwise() |>
  mutate(
    odds_ratio = unname(fisher.test(matrix(c(temperate_with, temperate_without,
                                             virulent_with, virulent_without),
                                           nrow = 2, byrow = TRUE))$estimate),
    fisher_p   = fisher.test(matrix(c(temperate_with, temperate_without,
                                      virulent_with, virulent_without),
                                    nrow = 2, byrow = TRUE))$p.value
  ) |>
  ungroup() |>
  mutate(fisher_p_adj_BH = p.adjust(fisher_p, method = "BH")) |>
  select(`function`, temperate_with, virulent_with, prop_temperate, prop_virulent,
         odds_ratio, fisher_p, fisher_p_adj_BH) |>
  arrange(fisher_p)

message(sprintf("=== Per-function Fisher's exact tests, temperate (n = %d) vs virulent (n = %d) vOTUs, BH over %d tests ===",
                n_temp_votus, n_vir_votus, nrow(lifestyle_function_fisher)))
print(as.data.frame(lifestyle_function_fisher |>
                      mutate(across(c(prop_temperate, prop_virulent, odds_ratio), ~ round(.x, 2)),
                             across(c(fisher_p, fisher_p_adj_BH), ~ signif(.x, 3)))),
      row.names = FALSE)
message("CAVEAT: this pools all 3 participants (lifestyle and function are fixed per-vOTU properties, not per-timepoint). The odds ratio is temperate:virulent (>1 = function more often carried by temperate vOTUs). Read alongside the uncalled vOTUs (about a third of the catalogue have no PhaTYP call).")

# -----------------------------------------------------------------------
# 4.6. AIM B FORMAL TEST: did AMG-carrying vOTU abundance actually change
#      around L1's antibiotic course, or does the trend line just look
#      that way?
# -----------------------------------------------------------------------
# Same pre- (timepoint <=4) vs post- (timepoint >=5) window split used
# throughout (Parts 2, 7, 8, 9). For L1 this is the genuine antibiotic
# window; for D3/G2 the identical split is only a structurally matched
# reference (no perturbation event for them), included for the same
# reason Part 8 included one - to judge whether L1's result is any bigger
# than the pre/post difference an unperturbed participant shows on its
# own. Unpaired Wilcoxon rank-sum test (samples pre and post aren't
# paired observations of the same unit).

amg_window_test <- function(part) {
  d <- amg_over_time |>
    filter(participant == part) |>
    mutate(window = if_else(timepoint <= 4, "pre", "post"))
  
  n_pre  <- sum(d$window == "pre")
  n_post <- sum(d$window == "post")
  
  if (n_pre < 2 || n_post < 2) {
    return(tibble(participant = part, n_pre = n_pre, n_post = n_post,
                  wilcoxon_p = NA_real_,
                  note = "too few AMG-positive samples in one window to test"))
  }
  
  wt <- suppressWarnings(wilcox.test(summed_abundance ~ window, data = d))
  tibble(participant = part, n_pre = n_pre, n_post = n_post,
         wilcoxon_p = wt$p.value, note = NA_character_)
}

amg_window_tests <- purrr::map_dfr(all_participants, amg_window_test)

message("=== Aim B: AMG-carrying vOTU abundance, pre- vs post-antibiotic-window Wilcoxon test (L1 = real event; D3/G2 = structurally matched reference only) ===")
print(amg_window_tests)
message("CAVEAT: n_pre/n_post here are only the samples where an AMG-carrying vOTU was actually detected at all (amg_over_time excludes zero-abundance samples), so these can be smaller than 4/5 and the Wilcoxon test has very limited power at this size - treat any p-value here as suggestive at best, and read it against the D3/G2 reference rows the same way Part 8's resilience check was read against reference_drift_check.")


# #############################################################################
# PART 5: CO-OCCURRENCE
# (originally aim_c_cooccurrence.R)
# Exploratory within-participant vOTU co-occurrence, based on CLR-
# transformed Spearman correlation.
#
# IMPORTANT SCOPE NOTE: each participant only has 8-10 timepoints. This
# is enough to compute pairwise correlations, but far too few for formal
# significance testing across ~200-1000 candidate pairs per participant
# (multiple-testing correction would leave little to nothing standing,
# not because the biology is absent but from lack of power). This
# analysis is deliberately EXPLORATORY/descriptive: edges are shown by
# correlation strength, not by a significance threshold. Frame it that
# way in the thesis - "hypothesis-generating," not confirmatory.
#
# Co-occurrence is computed WITHIN participant only (never pooling across
# people), matching the study's within-individual temporal focus.
# #############################################################################

# -----------------------------------------------------------------------
# 5.1. CLR transform helper
# -----------------------------------------------------------------------
# Addresses compositionality (abundance data is relative, not absolute -
# same closure problem discussed for normalization). Zeros are replaced
# with a small pseudocount (half the smallest observed non-zero value)
# rather than using a full zero-replacement package like zCompositions,
# to avoid another dependency for what is an exploratory analysis.
#
# Reused by Part 7 for phage:host pairing framing/caveats (though Part 7's
# actual correlation is computed on raw values, not CLR - see Section 7.5
# for why that distinction is appropriate there).

clr_transform_matrix <- function(mat) {
  # mat: vOTUs (rows) x samples (columns), numeric, >= 0
  pseudocount <- min(mat[mat > 0], na.rm = TRUE) / 2
  mat[mat == 0] <- pseudocount
  log_mat <- log(mat)
  sweep(log_mat, 2, colMeans(log_mat), FUN = "-")  # centre each sample (column) on its own geometric mean
}

# -----------------------------------------------------------------------
# 5.2. Candidate pairs - vOTU pairs with >=3 shared non-zero samples,
#      within participant (computed earlier during feasibility checking;
#      recomputed here so this section is self-contained)
# -----------------------------------------------------------------------

votu_sample_presence <- viral_long |>
  filter(abundance > 0) |>
  select(contig, sample_id, participant)

votu_pairs_strong <- votu_sample_presence |>
  inner_join(votu_sample_presence, by = c("sample_id", "participant"),
             relationship = "many-to-many") |>
  filter(contig.x < contig.y) |>
  count(contig.x, contig.y, participant, name = "shared_samples") |>
  filter(shared_samples >= 3)

message(sprintf("Candidate pairs (>=3 shared samples), by participant:"))
print(votu_pairs_strong |> count(participant))

# -----------------------------------------------------------------------
# 5.3. Per-participant CLR + correlation
# -----------------------------------------------------------------------

# Same-parent pairs: CheckV prophage regions cut from one assembled contig
# carry the suffix _1, _2, ... (e.g. L1.6__contig_2129_1 / _2). Two regions
# of one parent contig sit on the same physical molecule, so their read
# coverage is coupled by construction (same cell/chromosome, same mapping
# artefacts) and a co-occurrence between them is not independent evidence.
# Such pairs are dropped from the co-occurrence edges and the windowed
# networks. Both regions stay in the catalogue as separate vOTUs.
parent_contig_of <- function(x) str_replace(x, "^(.*__contig_[0-9]+)_[0-9]+$", "\\1")
is_same_parent <- function(a, b) parent_contig_of(a) == parent_contig_of(b)

compute_participant_cooccurrence <- function(part) {
  samples_this_participant <- viral_long |>
    filter(participant == part) |>
    distinct(sample_id) |>
    pull(sample_id)
  
  mat <- abundance_wide |>
    select(contig, all_of(samples_this_participant)) |>
    column_to_rownames("contig") |>
    as.matrix()
  
  # Drop vOTUs never detected in this participant's samples
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  
  clr_mat <- clr_transform_matrix(mat)
  cor_mat <- cor(t(clr_mat), method = "spearman")
  
  cor_mat |>
    as_tibble(rownames = "contig.x") |>
    pivot_longer(-contig.x, names_to = "contig.y", values_to = "correlation") |>
    filter(contig.x < contig.y) |>
    mutate(participant = part)
}

all_participants <- unique(viral_long$participant)
cooccurrence_all <- purrr::map_dfr(all_participants, compute_participant_cooccurrence)

# Restrict to the candidate pairs with real shared-sample support
cooccurrence_filtered <- cooccurrence_all |>
  inner_join(votu_pairs_strong, by = c("contig.x", "contig.y", "participant"))

same_parent_removed <- cooccurrence_filtered |>
  filter(is_same_parent(contig.x, contig.y))
message(sprintf("=== Same-parent prophage-region pairs removed from the co-occurrence edges: %d ===",
                nrow(same_parent_removed)))
if (nrow(same_parent_removed) > 0) print(same_parent_removed |> select(participant, contig.x, contig.y, correlation, shared_samples))
cooccurrence_filtered <- cooccurrence_filtered |>
  filter(!is_same_parent(contig.x, contig.y))

message(sprintf("Final edge count (>=3 shared samples, CLR Spearman correlation): %d",
                nrow(cooccurrence_filtered)))

# Clean vOTU display labels for the heatmap below (see Part 1, Section
# 1.1f) - contig.x/contig.y are KEPT as-is for the diagnostic tables
# further down (still the most useful form for tracing a specific pair
# back to its real contig), only the plot uses votu_display.x/.y.
cooccurrence_filtered <- cooccurrence_filtered |>
  left_join(votu_display_lookup |> select(contig, votu_display.x = votu_display),
            by = c("contig.x" = "contig")) |>
  left_join(votu_display_lookup |> select(contig, votu_display.y = votu_display),
            by = c("contig.y" = "contig"))

# -----------------------------------------------------------------------
# 5.3b. SMALL-N RELIABILITY TAG (bookmarked 2026-09-28, addressed here)
# -----------------------------------------------------------------------
# At the minimum floor (exactly 3 shared samples), a Spearman correlation
# has only 3! = 6 possible rank orderings, so it lands on one of a
# handful of discrete values and hits exactly +-1 whenever the ranks
# happen to agree/disagree completely - which is common at n=3 by chance
# alone, not because the two vOTUs are actually tightly coupled. This
# doesn't get worse gradually; it's a step change right at n=3, so
# rather than raising the floor for everyone (which would discard real,
# usable n=4+ pairs) or silently keeping n=3 pairs unflagged, every edge
# is tagged with a reliability tier below, and BOTH the full set and an
# n>=4-only version are shown - compare the two before treating any
# single n=3 "strong" correlation as a real finding.

cooccurrence_filtered <- cooccurrence_filtered |>
  mutate(reliability = case_when(
    shared_samples == 3 ~ "low (n=3 - correlation is discretized, can hit +-1 by chance)",
    shared_samples == 4 ~ "moderate (n=4)",
    TRUE                ~ "higher (n>=5)"
  ))

message("=== Co-occurrence edges by shared-sample reliability tier ===")
print(cooccurrence_filtered |> count(reliability))

# -----------------------------------------------------------------------
# 5.4. Visualization: correlation heatmap per participant
# -----------------------------------------------------------------------
# A correlation-magnitude heatmap rather than a formal network graph -
# avoids requiring igraph/ggraph, and communicates the same information
# (which pairs co-vary, and how strongly) without implying statistical
# confirmation this sample size can't support. n=3 (low-reliability)
# tiles are drawn at reduced opacity so the eye isn't drawn to them the
# same way as the more trustworthy n>=4 edges, without deleting them
# from the picture entirely.

for (part in all_participants) {
  part_data <- cooccurrence_filtered |> filter(participant == part)
  
  if (nrow(part_data) == 0) {
    message(sprintf("No qualifying pairs for %s - skipping plot", part))
    next
  }
  
  n_labels <- length(union(part_data$votu_display.x, part_data$votu_display.y))
  
  heatmap_dim <- max(8, n_labels * 0.16)
  
  p <- ggplot(part_data, aes(x = votu_display.x, y = votu_display.y, fill = correlation, alpha = shared_samples == 3)) +
    geom_tile() +
    scale_fill_gradient2(low = COL_DIV_LOW, mid = COL_DIV_MID, high = COL_DIV_HIGH,
                         midpoint = 0, limits = c(-1, 1)) +
    scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 1), guide = "none") +
    labs(title = sprintf("%s: within-participant vOTU co-occurrence (exploratory)", part),
         subtitle = "CLR-transformed Spearman correlation, pairs with >=3 shared samples - not significance-tested.\nFaded tiles = n=3 shared samples only (low reliability, see 5.3b). vOTU_### labels - see votu_display_lookup for the real contig IDs.",
         x = NULL, y = NULL, fill = "Correlation") +
    theme_thesis() +
    heatmap_axis_theme(n_labels) +
    theme(legend.position = "right")
  
  # width/height scale with how many vOTUs are on each axis - the image
  # saved to FIGURES_DIR (by style_figure(), see Part 2) is the one to
  # actually use in the thesis; the Plots-pane preview will still look
  # cramped for any participant with a lot of qualifying vOTUs (L1
  # especially).
  print(style_figure(p, width = heatmap_dim, height = heatmap_dim))
}

# Strongest pairs - shown two ways so the n=3 sensitivity is visible
# rather than hidden in one combined ranking.
message("=== Strongest pairs, ALL edges (n>=3 shared samples) - includes low-reliability n=3 pairs ===")
cooccurrence_filtered |>
  arrange(desc(abs(correlation))) |>
  slice_head(n = 15) |>
  select(participant, contig.x, contig.y, correlation, shared_samples, reliability) |>
  print(n = Inf)

message("=== Strongest pairs, RESTRICTED to n>=4 shared samples (the more trustworthy subset) ===")
cooccurrence_filtered |>
  filter(shared_samples >= 4) |>
  arrange(desc(abs(correlation))) |>
  slice_head(n = 15) |>
  select(participant, contig.x, contig.y, correlation, shared_samples, reliability) |>
  print(n = Inf)
message("Compare the two tables above: any pair whose correlation is exactly +-1.000 AND only appears in the n>=3 (not n>=4) list should be treated as a likely small-n artefact, not a real co-occurrence finding.")


# #############################################################################
# PART 6: VOLATILE vOTU FUNCTIONAL CONTEXT
# What do the "most variable" vOTUs (candidate_votus, Part 3 Section 3.7)
# actually carry, functionally? candidate_votus only carries a taxon label
# (Part 3) - this joins in real gene-level function calls from Phold
# (Part 4), so the volatility panels aren't just "this changed a lot" with
# no biological hook for why that might matter.
# #############################################################################

# -----------------------------------------------------------------------
# 6.1. Join candidate_votus -> votu_id -> Phold gene functions
# -----------------------------------------------------------------------
# candidate_votus is keyed by `contig` (the representative contig used for
# read-mapping abundance), not `votu_id` directly. viral_long already
# carries both (joined from confidence_tiers in Part 1), so it's the
# lookup table rather than rejoining confidence_tiers again (see the
# rejoin-collision note in Part 4, Section 4.3).

votu_id_lookup <- viral_long |> distinct(contig, votu_id)

candidate_votu_functions <- candidate_votus |>
  left_join(votu_id_lookup, by = "contig") |>
  left_join(taxonomy |> select(contig, taxon_label), by = "contig") |>
  left_join(
    phold_votu |> distinct(votu_id, `function`) |> filter(`function` != "unknown function"),
    by = "votu_id",
    relationship = "many-to-many"  # a vOTU can carry more than one function category
  )

message("=== Functional context for the 'most variable' vOTUs (candidate_votus) ===")
candidate_votu_functions |>
  group_by(participant, contig, taxon_label, cv, n_detected) |>
  summarise(
    functions = if (all(is.na(`function`))) "none annotated / no genes called"
    else paste(sort(unique(`function`)), collapse = "; "),
    .groups = "drop"
  ) |>
  arrange(participant, desc(cv)) |>
  print(n = Inf)

# NOTE: if you want this baked into the Part 3 Section 3.7 line-chart facet
# labels themselves (alongside taxon_label), that's a small follow-up
# change to votu_line_data's facet_label construction - left as a table
# here for now since several vOTUs carry multiple functions and cramming
# all of them into a facet strip gets crowded fast.


# #############################################################################
# PART 7: PHAGE:HOST RELATIONSHIP DYNAMICS
# PhaBOX2's predicted host (viral_long$Host) is SPECIES-level, formatted as
# "species:<binomial name>" (confirmed 2026-09-25: 9 distinct real host
# species predicted across the catalogue, e.g. "species:Faecalibacterium
# prausnitzii" - plus one generic, uninformative "species:bacterium" call,
# which is excluded below since it can't be matched to anything specific).
# This is a PREDICTED host (PhaBOX2/CHERRY-style host prediction), not a
# confirmed infection - frame findings here as candidate relationships, not
# established phage-host interactions.
# #############################################################################

# -----------------------------------------------------------------------
# 7.1. Extract predicted host species name and build candidate pair list
# -----------------------------------------------------------------------

viral_long <- viral_long |>
  mutate(
    host_species = if_else(
      Host != "-" & Host != "species:bacterium",
      str_remove(Host, "^species:"),
      NA_character_
    )
  )

# BUG FIX (found 2026-09-25, first run): viral_long has a zero-filled row
# for every contig in every sample (from the wide abundance matrix), so
# building candidate pairs straight off viral_long paired every
# host-predicted vOTU with ALL THREE participants, including participants
# where that vOTU was never actually detected. A phage abundance vector
# that's all zeros has no variance, so cor() correctly returned NA for
# those pairs downstream - that's why the first run showed 6 of 8 pairs
# as NA. Fixed by restricting to (contig, participant) combos where the
# vOTU was genuinely detected via mapping, using the same
# "detected >= N samples" floor as Part 5's co-occurrence pairs (>=3
# shared samples), so surviving pairs have enough real variance for a
# meaningful correlation, not just >0 detections.
phage_detected_by_participant <- viral_long |>
  filter(abundance > 0) |>
  count(contig, participant, name = "n_detected") |>
  filter(n_detected >= 3)

phage_host_candidates_raw <- viral_long |>
  filter(!is.na(host_species)) |>
  distinct(contig, participant, host_species)

phage_host_candidates <- phage_host_candidates_raw |>
  semi_join(phage_detected_by_participant, by = c("contig", "participant"))

message(sprintf(
  "=== Predicted phage-host candidate pairs: %d of %d survive the >=3-detections floor (dropped pairs were never actually observed in that participant) ===",
  nrow(phage_host_candidates), nrow(phage_host_candidates_raw)
))
print(phage_host_candidates |> count(host_species, sort = TRUE))

# -----------------------------------------------------------------------
# 7.2. Diagnostic: was the predicted host species actually detected via
#      Bracken in that participant's own samples? A predicted host that's
#      never observed in this dataset can't have a real abundance
#      relationship tracked, whatever PhaBOX2 predicted.
# -----------------------------------------------------------------------

bacterial_species_present <- bacterial_long |>
  distinct(bacterial_taxon, participant)

phage_host_pairs <- phage_host_candidates |>
  inner_join(bacterial_species_present,
             by = c("host_species" = "bacterial_taxon", "participant"))

message(sprintf(
  "=== Phage-host pairs with the predicted host actually detected: %d of %d candidate pairs ===",
  nrow(phage_host_pairs), nrow(phage_host_candidates)
))
print(phage_host_pairs)

if (nrow(phage_host_pairs) == 0) {
  message("No phage-host pairs have a detected host in this dataset - Sections 7.3-7.6 will produce no output. Worth re-checking host_species spelling/format against bacterial_long$bacterial_taxon if this is unexpected.")
}

# 7.2b. For predicted hosts that did NOT confirm - check whether that's a
# real absence or just a naming/strain-suffix mismatch between PhaBOX2's
# species string and Bracken's taxon name (e.g. "sp. <strain>" identifiers
# often don't round-trip exactly between databases/classifiers). Checked
# 2026-09-25: of 6 candidate pairs, only 2 confirmed - worth ruling this
# out before treating the other 4 as genuinely absent.
unmatched_hosts <- phage_host_candidates |>
  anti_join(bacterial_species_present, by = c("host_species" = "bacterial_taxon", "participant")) |>
  distinct(host_species) |>
  pull(host_species)

if (length(unmatched_hosts) > 0) {
  message("=== Unmatched predicted hosts - checking for a genus-level near-match in bacterial_long before concluding real absence ===")
  for (h in unmatched_hosts) {
    genus_word <- word(h, 1)
    near_matches <- bacterial_long |>
      filter(str_detect(bacterial_taxon, fixed(genus_word))) |>
      distinct(bacterial_taxon) |>
      pull(bacterial_taxon)
    message(sprintf("  '%s' - %d near-match(es) at genus level: %s",
                    h, length(near_matches),
                    if (length(near_matches) == 0) "none (genuinely absent from bacterial_long)"
                    else paste(near_matches, collapse = "; ")))
  }
}


# -----------------------------------------------------------------------
# 7.3. Paired abundance trajectories - phage vOTU vs predicted host,
#      min-max scaled within each pair so two completely different
#      abundance scales can be compared by SHAPE on the same axis, not
#      magnitude.
# -----------------------------------------------------------------------
# NORMALIZATION FIX (2026-09-29): previously used raw CoverM abundance
# for the phage and raw Bracken new_est_reads for the host - both are
# confounded by sequencing depth (a sample with more total reads gets
# higher raw values for everything in it, independent of true biological
# abundance), which is a different problem from the compositional-closure
# issue CLR addresses elsewhere in this file (Parts 5/8/9). Fixed by
# using GMPR-normalized viral abundance (abundance / size_factor, Part 1
# Section 1.3) for the phage side, and each sample's Bracken RELATIVE
# abundance fraction (that taxon's reads / that sample's total reads,
# post-contamination-filtering) for the host side - both now on a
# per-sample-depth-corrected basis before the shape comparison and
# correlation below.

bacterial_long <- bacterial_long |>
  group_by(sample_id) |>
  mutate(bracken_rel_abundance = new_est_reads / sum(new_est_reads)) |>
  ungroup()

minmax_scale <- function(x) {
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(rep(0, length(x)))
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

phage_abund <- viral_long |>
  mutate(phage_abund = abundance / size_factor) |>
  select(contig, participant, sample_id, timepoint, phage_abund)

host_abund <- bacterial_long |>
  select(bacterial_taxon, participant, sample_id, timepoint, host_abund = bracken_rel_abundance)

pair_trajectories <- phage_host_pairs |>
  left_join(phage_abund, by = c("contig", "participant")) |>
  left_join(host_abund,
            by = c("host_species" = "bacterial_taxon", "participant", "sample_id", "timepoint")) |>
  mutate(host_abund = replace_na(host_abund, 0)) |>
  left_join(sample_dates, by = "sample_id")

pair_plot_data <- pair_trajectories |>
  left_join(votu_display_lookup |> select(contig, votu_display), by = "contig") |>
  mutate(facet_label = paste0(votu_display, " <-> ", host_species)) |>
  group_by(facet_label) |>
  mutate(
    phage_scaled = minmax_scale(phage_abund),
    host_scaled  = minmax_scale(host_abund)
  ) |>
  ungroup() |>
  select(participant, contig, host_species, facet_label, sample_id, timepoint, collection_date,
         phage_scaled, host_scaled) |>
  pivot_longer(cols = c(phage_scaled, host_scaled), names_to = "series", values_to = "scaled_abundance") |>
  mutate(series = if_else(series == "phage_scaled", "Phage (vOTU)", "Predicted host"))

for (p in unique(pair_plot_data$participant)) {
  
  plot_data <- pair_plot_data |> filter(participant == p)
  if (nrow(plot_data) == 0) next
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = scaled_abundance, colour = series)) +
    geom_line() +
    geom_point(size = 1.5) +
    facet_wrap(~ facet_label, scales = "free_x") +
    scale_colour_manual(values = c("Phage (vOTU)" = COL_VIRAL, "Predicted host" = COL_BACTERIAL)) +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Phage:predicted-host abundance trajectories - ", p),
         subtitle = "GMPR-normalized phage abundance vs Bracken relative-abundance host fraction, min-max scaled within each pair for shape comparison",
         x = "Collection date", y = "Scaled abundance (0-1)", colour = NULL) +
    theme_thesis() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 7.4. Relationship-dynamics statistic: does the phage and its predicted
#      host co-vary over time?
# -----------------------------------------------------------------------
# Spearman correlation on raw values (NOT CLR-transformed, unlike Part 5's
# phage-phage screen). Part 5's CLR correction addresses compositionality
# when screening many candidate pairs at once, where one feature's
# apparent rise can be an artefact of another feature's fall within the
# same closed composition. Here we're asking a single targeted question
# about two already biologically-linked, independently-measured entities
# (phage abundance from CoverM, host abundance from Bracken - different
# compositions entirely), so that particular correction doesn't apply the
# same way; plain Spearman on the paired trajectories is the more
# interpretable choice for this comparison.

pair_correlations <- pair_trajectories |>
  mutate(facet_label = paste0(contig, " <-> ", host_species)) |>
  group_by(participant, contig, host_species, facet_label) |>
  summarise(
    n_timepoints = n(),
    correlation  = suppressWarnings(cor(phage_abund, host_abund, method = "spearman")),
    .groups = "drop"
  )

message("=== Phage:predicted-host correlation (exploratory - very small n per pair) ===")
print(pair_correlations, n = Inf)

# -----------------------------------------------------------------------
# 7.5. L1 specifically: does the phage-host relationship shift around the
#      antibiotic perturbation (timepoints 1-4 vs 5+)?
# -----------------------------------------------------------------------
# Exploratory/descriptive only - each window has very few timepoints (L1
# has 9 samples total, so 4 pre / 5 post at most), so treat any apparent
# shift as a candidate finding worth highlighting, not a tested claim.
#
# ZERO-VARIANCE FIX (2026-09-29): the full-series correlation (Section
# 7.4) already guards against a pair with no variance across all 9
# timepoints, but that guarantee doesn't carry over to a smaller 4- or
# 5-point sub-window specifically - a pair can vary across the full
# series while still being constant WITHIN one window (e.g.
# L1.6__contig_1404 / L. delbrueckii was constant across just the
# pre-antibiotic samples), which cor() silently returns as NA for with
# no explanation. Checked explicitly below (phage_var_zero/host_var_zero)
# so a NA in this table is now labelled with WHY, not left to guesswork.

l1_pair_window <- pair_trajectories |>
  filter(participant == "L1") |>
  mutate(window = if_else(timepoint <= 4, "pre-antibiotic", "post-antibiotic"))

if (nrow(l1_pair_window) > 0) {
  message("=== L1 phage:predicted-host correlation, pre vs post antibiotic window (exploratory) ===")
  l1_pair_window |>
    group_by(contig, host_species, window) |>
    summarise(
      n_timepoints   = n(),
      phage_var_zero = length(unique(phage_abund)) <= 1,
      host_var_zero  = length(unique(host_abund)) <= 1,
      correlation    = if (phage_var_zero[1] || host_var_zero[1]) {
        NA_real_
      } else {
        suppressWarnings(cor(phage_abund, host_abund, method = "spearman"))
      },
      na_reason = case_when(
        phage_var_zero[1] && host_var_zero[1] ~ "both phage and host abundance constant within this window",
        phage_var_zero[1] ~ "phage abundance constant within this window",
        host_var_zero[1]  ~ "host abundance constant within this window",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) |>
    arrange(contig, host_species, window) |>
    print(n = Inf)
} else {
  message("No phage-host pairs detected for L1 - no pre/post antibiotic comparison to show.")
}

# #############################################################################
# PART 8: COMMUNITY RESILIENCE
# Does community composition move back toward its pre-antibiotic state after
# L1's perturbation, or settle into a new one? Previously the antibiotic
# window was only ever marked visually (antibiotic_marker()) on trend lines -
# this is the actual metric behind that marker: Aitchison distance
# (Euclidean distance on CLR-transformed abundances) from a baseline
# centroid, tracked per timepoint. Reuses Part 5's clr_transform_matrix()
# for methodological consistency, rather than introducing a second distance
# metric (e.g. Bray-Curtis) with different assumptions elsewhere in the
# analysis.
# #############################################################################

# -----------------------------------------------------------------------
# 8.1. Per-participant CLR matrices, viral and bacterial
# -----------------------------------------------------------------------
# Viral: same abundance_wide subsetting Part 5 uses for co-occurrence.
# Bacterial: bacterial_long has no wide-matrix equivalent yet, so one is
# built here per participant (taxa x samples, zero-filled), then CLR
# transformed with Part 5's helper - same treatment, different domain.

get_participant_clr <- function(part, domain = c("viral", "bacterial")) {
  domain <- match.arg(domain)
  
  if (domain == "viral") {
    samples_this_participant <- viral_long |>
      filter(participant == part) |>
      distinct(sample_id) |>
      pull(sample_id)
    
    mat <- abundance_wide |>
      select(contig, all_of(samples_this_participant)) |>
      column_to_rownames("contig") |>
      as.matrix()
  } else {
    mat <- bacterial_long |>
      filter(participant == part) |>
      select(bacterial_taxon, sample_id, new_est_reads) |>
      pivot_wider(names_from = sample_id, values_from = new_est_reads,
                  values_fill = 0) |>
      column_to_rownames("bacterial_taxon") |>
      as.matrix()
  }
  
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  clr_transform_matrix(mat)
}

# -----------------------------------------------------------------------
# 8.2. Distance-from-baseline: Aitchison distance between each timepoint's
#      CLR profile and the centroid of a defined "baseline" window.
# -----------------------------------------------------------------------
# For L1, baseline = the 4 pre-antibiotic timepoints - a genuinely
# meaningful pre-perturbation state. For D3/G2, the SAME timepoint<=4 split
# is applied purely as a STRUCTURALLY MATCHED reference window - there is
# no perturbation event for them, so their trajectory carries no resilience
# meaning of its own. It exists only to show what ordinary background
# compositional drift looks like on this same metric, so L1's
# post-antibiotic distances can be judged against a real magnitude
# reference instead of in isolation.

distance_from_baseline <- function(clr_mat, baseline_samples) {
  baseline_samples <- intersect(baseline_samples, colnames(clr_mat))
  if (length(baseline_samples) == 0) return(NULL)
  centroid <- rowMeans(clr_mat[, baseline_samples, drop = FALSE])
  dists <- apply(clr_mat, 2, function(s) sqrt(sum((s - centroid)^2)))
  tibble(sample_id = names(dists), aitchison_distance = dists)
}

compute_resilience_trajectory <- function(part, domain) {
  clr_mat <- get_participant_clr(part, domain)
  
  samples_this_participant <- viral_long |>
    filter(participant == part) |>
    distinct(sample_id, timepoint, collection_date) |>
    arrange(timepoint)
  
  baseline_samples <- samples_this_participant |>
    filter(timepoint <= 4) |>
    pull(sample_id)
  
  distance_from_baseline(clr_mat, baseline_samples) |>
    left_join(samples_this_participant, by = "sample_id") |>
    mutate(participant = part, domain = domain)
}

resilience_trajectory <- purrr::map_dfr(all_participants, function(p) {
  bind_rows(
    compute_resilience_trajectory(p, "viral"),
    compute_resilience_trajectory(p, "bacterial")
  )
})

message("=== Community distance-from-baseline (Aitchison/CLR), all participants ===")
print(resilience_trajectory |> arrange(domain, participant, timepoint), n = Inf)

# -----------------------------------------------------------------------
# 8.3. Visualization: distance-from-baseline over time, per participant,
#      viral vs bacterial, with L1's antibiotic window marked.
# -----------------------------------------------------------------------

p_resilience_trajectory <- ggplot(resilience_trajectory,
                                  aes(x = collection_date, y = aitchison_distance, colour = domain)) +
  geom_line() +
  geom_point(size = 2) +
  antibiotic_marker_date() +
  facet_wrap(~ participant, scales = "free_x") +
  scale_colour_manual(values = c("viral" = COL_VIRAL, "bacterial" = COL_BACTERIAL)) +
  participant_date_scales() +
  labs(title = "Community distance from baseline over time (Aitchison/CLR distance)",
       subtitle = "Baseline = centroid of each participant's first 4 timepoints. For L1 this is the genuine pre-antibiotic state; for D3/G2 it's a structurally matched reference window only (no perturbation event), shown to gauge normal background drift",
       x = "Collection date", y = "Aitchison distance from baseline centroid", colour = NULL) +
  theme_thesis() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

print(style_figure(p_resilience_trajectory))

# -----------------------------------------------------------------------
# 8.4. Headline resilience check: does distance from baseline decline from
#      its PEAK back toward baseline, by the last available timepoint?
# -----------------------------------------------------------------------
# REVISED 2026-09-26 (first run used first-vs-last, not peak-vs-last): the
# first-vs-last comparison silently assumes the perturbation's maximal
# compositional effect lands on the very first post-window sample. On the
# real data it doesn't - both of L1's domains peak at timepoint 6 (one
# sample after the antibiotic window), not timepoint 5 - so first-vs-last
# understated the viral decline entirely (it looked flat) by comparing two
# points that both happened to sit below the actual peak. Peak-vs-last is
# the standard ecological framing anyway: peak = resistance/impact
# magnitude, decline-from-peak = resilience.

l1_resilience_check <- resilience_trajectory |>
  filter(participant == "L1", timepoint >= 5) |>
  group_by(domain) |>
  summarise(
    peak_distance         = max(aitchison_distance),
    peak_timepoint        = timepoint[which.max(aitchison_distance)],
    last_distance         = aitchison_distance[which.max(timepoint)],
    last_timepoint        = max(timepoint),
    change_from_peak      = last_distance - peak_distance,
    pct_change_from_peak  = round(100 * change_from_peak / peak_distance, 1),
    direction             = if_else(change_from_peak < 0,
                                    "declining from peak toward baseline (resilience)",
                                    "still at or above peak (no recovery within this window)"),
    .groups = "drop"
  )

message("=== L1 resilience check: distance from baseline, peak vs last post-antibiotic timepoint ===")
print(l1_resilience_check)

# For context: does D3/G2 (no perturbation) show a comparable peak-then-
# decline shape over an equivalent span? If so, a decline alone isn't
# strong evidence the antibiotic specifically drove L1's recovery pattern,
# rather than ordinary compositional random walk / regression-to-the-mean
# behaviour that shows up in unperturbed samples too.
reference_drift_check <- resilience_trajectory |>
  filter(participant != "L1", timepoint >= 5) |>
  group_by(participant, domain) |>
  summarise(
    peak_distance    = max(aitchison_distance),
    peak_timepoint   = timepoint[which.max(aitchison_distance)],
    last_distance    = aitchison_distance[which.max(timepoint)],
    last_timepoint   = max(timepoint),
    change_from_peak = last_distance - peak_distance,
    .groups = "drop"
  )

message("=== Reference: same peak-vs-last comparison for D3/G2 (no perturbation - background drift magnitude) ===")
print(reference_drift_check)
message("NOTE: if a reference participant shows a decline-from-peak of comparable size to L1's, that specific participant's data can't tell resilience apart from ordinary drift. The disruption itself (l1_noise_vs_shift, pre vs post mean distance) is the more robust finding here - the peak-to-last decline is suggestive, not conclusive, given only one perturbed individual and ~5 post-antibiotic timepoints.")

# -----------------------------------------------------------------------
# 8.4b. L1: average distance from baseline, pre-window (self-variability /
#       "noise floor") vs post-window (the shift itself) - gives a sense of
#       whether the post-antibiotic distances are meaningfully bigger than
#       L1's own normal sample-to-sample wobble, not just bigger than zero.
# -----------------------------------------------------------------------

l1_noise_vs_shift <- resilience_trajectory |>
  filter(participant == "L1") |>
  mutate(window = if_else(timepoint <= 4, "pre-antibiotic (baseline noise)", "post-antibiotic")) |>
  group_by(domain, window) |>
  summarise(mean_distance = mean(aitchison_distance), .groups = "drop")

message("=== L1: average distance from baseline, pre (self-variability) vs post (shift) ===")
print(l1_noise_vs_shift)

# #############################################################################
# PART 9: NETWORK STABILITY
# Does L1's vOTU co-occurrence network (Part 5) hold its structure across the
# antibiotic perturbation, or reorganise? Scope is L1 only - the only
# participant with an actual perturbation event to test stability against.
# Reuses Part 5's clr_transform_matrix() and CLR-Spearman approach, applied
# separately to the pre-antibiotic window (timepoints 1-4) and the
# post-antibiotic window (timepoints 5+), then the two windows' networks are
# compared: how many nodes/edges, how dense, and how much of the edge set is
# actually the SAME set of pairs (Jaccard overlap) vs a reorganisation.
# #############################################################################

# -----------------------------------------------------------------------
# 9.1. Windowed network builder
# -----------------------------------------------------------------------
# min_shared is deliberately lower than Part 5's >=3-shared-sample threshold,
# which was applied across a participant's full 8-10 timepoints. A single
# window here only has 4-5 samples, so a fixed >=3 would be disproportionately
# strict relative to the window size. min_shared = 2 keeps roughly the same
# PROPORTIONAL stringency (about half the samples in the window actually
# sharing detection) rather than being a leniency shortcut.
#
# "Eligible" nodes (for density) are vOTUs detected in >=min_shared samples
# WITHIN that window - not just vOTUs that happen to appear in an edge - so
# density isn't circularly defined by the edges themselves.
#
# SMALL-N CAVEAT (corrected 2026-10-14): the edge FILTER is "shared_samples
# >= min_shared", but the Spearman correlation itself is computed across ALL
# samples in the window (CLR values, zeros replaced), not only the shared
# ones, so a min_shared = 2 edge is NOT forced to be +-1 (the earlier comment
# here said it was; that was wrong). With only 4-5 samples per window,
# though, each correlation can take only a handful of values and is
# unstable, so edge counts, density and mean |r| remain inflated/noisy.
# Section 9.6 rebuilds both windows at min_shared = 3 as a sensitivity check
# and 9.6b restricts the edge-set comparison to vOTUs eligible in both
# windows and adds the D3/G2 reference.

build_window_network <- function(sample_ids, min_shared = 2) {
  mat <- abundance_wide |>
    select(contig, all_of(sample_ids)) |>
    column_to_rownames("contig") |>
    as.matrix()
  
  mat <- mat[rowSums(mat > 0) >= min_shared, , drop = FALSE]
  
  if (nrow(mat) < 2) {
    return(list(edges = tibble(contig.x = character(), contig.y = character(),
                               correlation = numeric(), shared_samples = integer()),
                eligible_votus = rownames(mat)))
  }
  
  clr_mat <- clr_transform_matrix(mat)
  cor_mat <- cor(t(clr_mat), method = "spearman")
  
  presence <- (mat > 0)
  shared_counts <- presence %*% t(presence)  # vOTU x vOTU shared-sample counts
  
  shared_long <- shared_counts |>
    as_tibble(rownames = "contig.x") |>
    pivot_longer(-contig.x, names_to = "contig.y", values_to = "shared_samples") |>
    filter(contig.x < contig.y)
  
  edges <- cor_mat |>
    as_tibble(rownames = "contig.x") |>
    pivot_longer(-contig.x, names_to = "contig.y", values_to = "correlation") |>
    filter(contig.x < contig.y) |>
    inner_join(shared_long, by = c("contig.x", "contig.y")) |>
    filter(shared_samples >= min_shared) |>
    filter(!is_same_parent(contig.x, contig.y))  # see parent_contig_of() in Part 5
  
  list(edges = edges, eligible_votus = rownames(mat))
}

# -----------------------------------------------------------------------
# 9.2. L1's pre- and post-antibiotic windows
# -----------------------------------------------------------------------

l1_samples_by_window <- viral_long |>
  filter(participant == "L1") |>
  distinct(sample_id, timepoint) |>
  arrange(timepoint) |>
  mutate(window = if_else(timepoint <= 4, "pre-antibiotic", "post-antibiotic"))

l1_pre_samples  <- l1_samples_by_window |> filter(window == "pre-antibiotic")  |> pull(sample_id)
l1_post_samples <- l1_samples_by_window |> filter(window == "post-antibiotic") |> pull(sample_id)

l1_network_pre  <- build_window_network(l1_pre_samples)
l1_network_post <- build_window_network(l1_post_samples)

# -----------------------------------------------------------------------
# 9.3. Network summary: nodes, edges, density, mean edge strength
# -----------------------------------------------------------------------

summarise_network <- function(net, window_label) {
  n_eligible         <- length(net$eligible_votus)
  n_edges            <- nrow(net$edges)
  max_possible_edges <- choose(n_eligible, 2)
  
  tibble(
    window               = window_label,
    n_eligible_votus     = n_eligible,
    n_edges              = n_edges,
    density              = if_else(max_possible_edges > 0, n_edges / max_possible_edges, NA_real_),
    mean_abs_correlation = if (n_edges > 0) mean(abs(net$edges$correlation)) else NA_real_
  )
}

network_summary <- bind_rows(
  summarise_network(l1_network_pre,  "pre-antibiotic (timepoints 1-4)"),
  summarise_network(l1_network_post, "post-antibiotic (timepoints 5+)")
)

message("=== L1 network summary: pre vs post antibiotic ===")
print(network_summary)

# -----------------------------------------------------------------------
# 9.4. Edge-set turnover: Jaccard overlap of pre vs post edge sets
# -----------------------------------------------------------------------
# Jaccard = |shared edges| / |union of edges|. 0 = the two windows share no
# edges at all (complete reorganisation); 1 = identical edge sets (fully
# stable). This is about which SPECIFIC pairs recur, not just whether the
# summary stats (9.3) happen to look similar - two windows can have the same
# edge count and density while being made of entirely different pairs.

edge_id <- function(edges_df) paste(edges_df$contig.x, edges_df$contig.y, sep = "__")

pre_edge_ids  <- edge_id(l1_network_pre$edges)
post_edge_ids <- edge_id(l1_network_post$edges)

shared_edges    <- intersect(pre_edge_ids, post_edge_ids)
pre_only_edges  <- setdiff(pre_edge_ids, post_edge_ids)
post_only_edges <- setdiff(post_edge_ids, pre_edge_ids)
union_edges     <- union(pre_edge_ids, post_edge_ids)

jaccard_index <- if (length(union_edges) > 0) length(shared_edges) / length(union_edges) else NA_real_

network_turnover <- tibble(
  n_pre_edges    = length(pre_edge_ids),
  n_post_edges   = length(post_edge_ids),
  n_shared_edges = length(shared_edges),
  n_pre_only     = length(pre_only_edges),
  n_post_only    = length(post_only_edges),
  jaccard_index  = jaccard_index
)

message("=== L1 edge-set turnover, pre vs post antibiotic (Jaccard index of shared edges) ===")
print(network_turnover)
message(sprintf(
  "Jaccard index = %s: %s (0 = complete edge-set turnover/full reorganisation, 1 = fully stable/identical edge sets). %d of %d union edges are shared between windows.",
  if (is.na(jaccard_index)) "NA" else sprintf("%.3f", jaccard_index),
  if (is.na(jaccard_index)) "not computable (no edges in one or both windows)"
  else if (jaccard_index >= 0.5) "more than half the network held together across the perturbation"
  else "the network largely reorganised across the perturbation",
  length(shared_edges), length(union_edges)
))

# -----------------------------------------------------------------------
# 9.5. Visualization: pre vs post network heatmaps, side by side
# -----------------------------------------------------------------------
# Same CLR-Spearman correlation-magnitude heatmap style as Part 5, so it
# reads as the same kind of evidence - exploratory, not significance-tested,
# given the very small per-window sample sizes (4-5 timepoints each).

network_edges_combined <- bind_rows(
  l1_network_pre$edges  |> mutate(window = "pre-antibiotic"),
  l1_network_post$edges |> mutate(window = "post-antibiotic")
) |>
  # Clean vOTU display labels (Part 1, Section 1.1f) - contig.x/contig.y
  # are kept alongside for anything downstream that still wants the real
  # contig ID (e.g. the Jaccard/edge-set diagnostics above).
  left_join(votu_display_lookup |> select(contig, votu_display.x = votu_display),
            by = c("contig.x" = "contig")) |>
  left_join(votu_display_lookup |> select(contig, votu_display.y = votu_display),
            by = c("contig.y" = "contig")) |>
  # facet_wrap orders panels by factor level (alphabetical otherwise, which
  # would put "post" on the left): force pre-antibiotic left, post right.
  mutate(window = factor(window, levels = c("pre-antibiotic", "post-antibiotic")))

if (nrow(network_edges_combined) > 0) {
  n_labels <- max(length(union(l1_network_pre$edges$contig.x,  l1_network_pre$edges$contig.y)),
                  length(union(l1_network_post$edges$contig.x, l1_network_post$edges$contig.y)))
  
  heatmap_dim <- max(8, n_labels * 0.16)
  
  p_network_stability <- ggplot(network_edges_combined,
                                aes(x = votu_display.x, y = votu_display.y, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = COL_DIV_LOW, mid = COL_DIV_MID, high = COL_DIV_HIGH,
                         midpoint = 0, limits = c(-1, 1)) +
    facet_wrap(~ window, scales = "free") +
    labs(title = "L1: vOTU co-occurrence network structure, pre vs post antibiotic (exploratory)",
         subtitle = sprintf("CLR-transformed Spearman correlation, pairs with >=2 shared samples within window - not significance-tested. Edge-set Jaccard overlap = %s. vOTU_### labels - see votu_display_lookup for the real contig IDs.",
                            if (is.na(jaccard_index)) "NA" else sprintf("%.3f", jaccard_index)),
         x = NULL, y = NULL, fill = "Correlation") +
    theme_thesis() +
    heatmap_axis_theme(n_labels) +
    theme(legend.position = "right")
  
  # width/height scale with how many vOTUs are on each axis - this is
  # the L1 co-occurrence figure flagged as unusable at its default size;
  # the two side-by-side (pre/post) panels need a wide canvas.
  print(style_figure(p_network_stability, width = heatmap_dim * 2, height = heatmap_dim))
} else {
  message("No edges in either window - skipping network stability plot.")
}

# -----------------------------------------------------------------------
# 9.6. SENSITIVITY CHECK: same pre/post comparison rebuilt at min_shared=3
#      (bookmarked 2026-09-28, addressed here)
# -----------------------------------------------------------------------
# Rebuilds both windows with the same Part-5-matching >=3-shared-sample
# floor used elsewhere in this file and reports the same summary stats next
# to the min_shared = 2 version above. Caution when reading the Jaccard:
# the number of eligible vOTUs differs between windows at this floor, which
# by itself lowers the overlap - see 9.6b for the shared-node version.

l1_network_pre_strict  <- build_window_network(l1_pre_samples,  min_shared = 3)
l1_network_post_strict <- build_window_network(l1_post_samples, min_shared = 3)

network_summary_strict <- bind_rows(
  summarise_network(l1_network_pre_strict,  "pre-antibiotic (timepoints 1-4)"),
  summarise_network(l1_network_post_strict, "post-antibiotic (timepoints 5+)")
)

pre_edge_ids_strict  <- edge_id(l1_network_pre_strict$edges)
post_edge_ids_strict <- edge_id(l1_network_post_strict$edges)
shared_edges_strict  <- intersect(pre_edge_ids_strict, post_edge_ids_strict)
union_edges_strict   <- union(pre_edge_ids_strict, post_edge_ids_strict)
jaccard_index_strict <- if (length(union_edges_strict) > 0) {
  length(shared_edges_strict) / length(union_edges_strict)
} else {
  NA_real_
}

message("=== SENSITIVITY CHECK: L1 network summary at min_shared=3 (vs min_shared=2 above) ===")
print(network_summary_strict)
message(sprintf(
  "min_shared=2 (Section 9.3/9.4): %d pre-edges, %d post-edges, Jaccard = %s.\nmin_shared=3 (this check):     %d pre-edges, %d post-edges, Jaccard = %s.",
  length(pre_edge_ids), length(post_edge_ids),
  if (is.na(jaccard_index)) "NA" else sprintf("%.3f", jaccard_index),
  length(pre_edge_ids_strict), length(post_edge_ids_strict),
  if (is.na(jaccard_index_strict)) "NA" else sprintf("%.3f", jaccard_index_strict)
))

# -----------------------------------------------------------------------
# 9.6b. EDGE-SET TURNOVER RESTRICTED TO SHARED NODES + D3/G2 REFERENCE
#       (added 2026-10-14)
# -----------------------------------------------------------------------
# Why: an edge can only be shared between windows if BOTH of its vOTUs are
# eligible (detected in >= min_shared samples) in BOTH windows. The plain
# Jaccard in 9.4/9.6 therefore falls partly because the two windows have
# different node sets (e.g. 15 vs 31 eligible vOTUs at min_shared = 3), not
# only because the same vOTUs re-wired. This block (a) restricts both edge
# sets to vOTUs eligible in both windows, (b) gives the Jaccard expected if
# edges were placed at random with each window's own density on that shared
# node set (a benchmark, ignores degree structure), (c) reports how many
# shared edges keep the same sign, and (d) repeats everything for D3 and G2
# with the same timepoint split (<=4 vs >=5; no perturbation), so L1's value
# can be read against ordinary window-to-window re-wiring.

window_turnover_for <- function(part, min_shared) {
  samp <- viral_long |>
    filter(participant == part) |>
    distinct(sample_id, timepoint)
  pre_ids  <- samp |> filter(timepoint <= 4) |> pull(sample_id)
  post_ids <- samp |> filter(timepoint >= 5) |> pull(sample_id)
  
  net_pre  <- build_window_network(pre_ids,  min_shared = min_shared)
  net_post <- build_window_network(post_ids, min_shared = min_shared)
  
  nodes_both <- intersect(net_pre$eligible_votus, net_post$eligible_votus)
  keep_edges <- function(e) e |> filter(contig.x %in% nodes_both, contig.y %in% nodes_both)
  e_pre  <- keep_edges(net_pre$edges)
  e_post <- keep_edges(net_post$edges)
  
  ids_pre  <- edge_id(net_pre$edges);  ids_post  <- edge_id(net_post$edges)
  sn_pre   <- edge_id(e_pre);          sn_post   <- edge_id(e_post)
  shared_sn <- intersect(sn_pre, sn_post)
  
  n_pairs <- choose(length(nodes_both), 2)
  d1 <- if (n_pairs > 0) length(sn_pre)  / n_pairs else NA_real_
  d2 <- if (n_pairs > 0) length(sn_post) / n_pairs else NA_real_
  
  sign_agree <- NA_real_
  if (length(shared_sn) > 0) {
    s_pre  <- setNames(sign(e_pre$correlation),  sn_pre)[shared_sn]
    s_post <- setNames(sign(e_post$correlation), sn_post)[shared_sn]
    sign_agree <- mean(s_pre == s_post)
  }
  
  tibble(
    participant = part, min_shared = min_shared,
    n_pre_samples = length(pre_ids), n_post_samples = length(post_ids),
    eligible_pre = length(net_pre$eligible_votus),
    eligible_post = length(net_post$eligible_votus),
    shared_nodes = length(nodes_both),
    jaccard_all_edges = if (length(union(ids_pre, ids_post)) > 0)
      length(intersect(ids_pre, ids_post)) / length(union(ids_pre, ids_post)) else NA_real_,
    edges_pre_sharednodes = length(sn_pre),
    edges_post_sharednodes = length(sn_post),
    edges_shared = length(shared_sn),
    jaccard_shared_nodes = if (length(union(sn_pre, sn_post)) > 0)
      length(shared_sn) / length(union(sn_pre, sn_post)) else NA_real_,
    jaccard_random_expectation = if (!is.na(d1) && (d1 + d2 - d1 * d2) > 0)
      (d1 * d2) / (d1 + d2 - d1 * d2) else NA_real_,
    shared_edges_same_sign = sign_agree
  )
}

network_turnover_reference <- tidyr::expand_grid(
  participant = c("L1", "D3", "G2"), min_shared = c(2, 3)
) |>
  purrr::pmap_dfr(function(participant, min_shared) window_turnover_for(participant, min_shared))

message("=== Edge-set turnover restricted to vOTUs eligible in BOTH windows, L1 (antibiotic) vs D3/G2 (same split, no event) ===")
print(as.data.frame(network_turnover_reference |>
                      mutate(across(where(is.double) & !c(min_shared), ~ round(.x, 3)))),
      row.names = FALSE)
message("READING THIS: jaccard_shared_nodes is the like-for-like overlap (same vOTUs available in both windows). jaccard_random_expectation is what independent random edges at the same densities would give on that node set. If L1's jaccard_shared_nodes is well below its random expectation AND clearly lower than D3/G2's, the re-wiring is specific to L1; if D3/G2 look the same, this is ordinary window-to-window turnover at 4-6 samples per window. Window sizes differ (D3 4 vs 6, G2 4 vs 4, L1 4 vs 5).")

# -----------------------------------------------------------------------
# 9.7. FORMAL TEST: did overall diversity/abundance magnitude actually
#      shift pre vs post antibiotic, or does the network's density/edge
#      summary (9.3) + Jaccard turnover (9.4) already cover this?
# -----------------------------------------------------------------------
# These are complementary, not redundant, questions. Density/mean-edge-
# strength/Jaccard (9.3-9.6) describe the network's RELATIONAL structure -
# whether the same vOTUs co-vary with each other before and after - which
# says nothing about whether the community's overall SIZE or DIVERSITY
# changed. Two windows could have near-identical network density while
# richness collapses, or vice versa. This section tests that separate,
# simpler question directly: same pre- (timepoint <=4) vs post- (timepoint
# >=5) window split used throughout, same unpaired Wilcoxon rank-sum
# approach as Part 4 Section 4.6's AMG test, and the same D3/G2
# structurally-matched-reference framing (no perturbation event for them -
# included only to show what an ordinary, unperturbed pre/post difference
# looks like on each metric, so L1's result can be read against that
# rather than in isolation). Six metrics per participant: viral and
# bacterial richness, Shannon diversity, and total abundance (viral:
# GMPR-normalized; bacterial: total Bracken reads per sample).
#
# No BH correction applied across these 18 tests (3 participants x 2
# domains x 3 metrics), matching the AMG test's precedent - this is an
# exploratory add-on comparison, not a new formal Aim-level claim the way
# Aim A's Kendall trend tests (Part 2, Section 2.6) are.

bacterial_total_abundance <- bacterial_long |>
  group_by(sample_id, participant, timepoint) |>
  summarise(total_reads = sum(new_est_reads), .groups = "drop")

window_test_input <- bind_rows(
  viral_diversity |>
    select(participant, timepoint, richness, shannon) |>
    pivot_longer(cols = c(richness, shannon), names_to = "metric", values_to = "value") |>
    mutate(domain = "Viral"),
  bacterial_diversity |>
    select(participant, timepoint, richness, shannon) |>
    pivot_longer(cols = c(richness, shannon), names_to = "metric", values_to = "value") |>
    mutate(domain = "Bacterial"),
  viral_total_abundance |>
    transmute(participant, timepoint, metric = "total_abundance",
              value = normalized_abundance, domain = "Viral"),
  bacterial_total_abundance |>
    transmute(participant, timepoint, metric = "total_abundance",
              value = total_reads, domain = "Bacterial")
)

window_wilcoxon_test <- function(d, ...) {
  d <- d |> mutate(window = if_else(timepoint <= 4, "pre", "post"))
  n_pre  <- sum(d$window == "pre")
  n_post <- sum(d$window == "post")
  
  if (n_pre < 2 || n_post < 2) {
    return(tibble(n_pre = n_pre, n_post = n_post, wilcoxon_p = NA_real_,
                  note = "too few samples in one window to test"))
  }
  if (length(unique(d$value)) <= 1) {
    return(tibble(n_pre = n_pre, n_post = n_post, wilcoxon_p = NA_real_,
                  note = "no variance in this metric across the full series - test undefined"))
  }
  
  wt <- suppressWarnings(wilcox.test(value ~ window, data = d))
  tibble(n_pre = n_pre, n_post = n_post, wilcoxon_p = wt$p.value, note = NA_character_)
}

diversity_abundance_window_tests <- window_test_input |>
  group_by(participant, domain, metric) |>
  group_modify(~ window_wilcoxon_test(.x)) |>
  ungroup() |>
  arrange(domain, metric, participant)

message("=== Viral/bacterial diversity + total abundance, pre- vs post-antibiotic-window Wilcoxon test (L1 = real event; D3/G2 = structurally matched reference only) ===")
print(diversity_abundance_window_tests, n = Inf)
message("CAVEAT: pre/post window sizes here are only 4-5 samples per participant (the full series, not a detected-only subset), so this test has very limited power at this size, same caveat as the AMG test in Part 4 Section 4.6 - treat any p-value here as suggestive at best, read against the D3/G2 reference rows rather than in isolation, and alongside Part 8's Aitchison-distance resilience trajectory (which tracks the SAME pre/post contrast as a continuous magnitude rather than a single test statistic).")


# #############################################################################
# PART 10: DESCRIPTIVE COMMUNITY COMPOSITION
# Purely descriptive figures (no hypothesis tests) answering "what specific
# bacteria and viruses were actually there, and in what abundance" - this
# is a mostly observational study, and Part 3's CV-selected "most/least
# variable" panels answer a different question (what changed the most/
# least), not "what dominates this person's gut community". This part
# selects by raw mean abundance instead of CV, and adds a single
# across-taxa heatmap view per participant/domain so the whole community
# picture can be read at a glance rather than assembled from 8 separate
# line charts. Reuses genus_reads_full (Part 3, Section 3.6),
# species_reads_full (Part 3, Section 3.6C), and votu_reads_full (Part 3,
# Section 3.7) - all already zero-filled across each participant's full
# sample set and already carry collection_date. Section 10.3B (species
# abundance heatmap, NEW 2026-09-30) mirrors 10.3's genus heatmap at
# species-level resolution.
# #############################################################################

# -----------------------------------------------------------------------
# 10.1. BACTERIAL - trajectories of the most ABUNDANT genera (by mean
#       reads), per participant - not CV-selected, so this shows the
#       community's dominant members regardless of how much they moved.
# -----------------------------------------------------------------------

top_abundant_genera <- genus_reads_full |>
  group_by(participant, bacterial_genus) |>
  summarise(mean_reads = mean(reads), .groups = "drop") |>
  group_by(participant) |>
  slice_max(mean_reads, n = 8) |>
  arrange(participant, desc(mean_reads))

message("=== Most abundant bacterial genera per participant (by mean reads, not CV) ===")
print(top_abundant_genera, n = Inf)

top_abundant_genus_line_data <- genus_reads_full |>
  semi_join(top_abundant_genera, by = c("participant", "bacterial_genus"))

for (p in unique(top_abundant_genus_line_data$participant)) {
  
  plot_data <- top_abundant_genus_line_data |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = reads)) +
    geom_line(colour = COL_BACTERIAL) +
    geom_point(size = 1.8, colour = COL_BACTERIAL) +
    facet_wrap(~ bacterial_genus, scales = "free_y") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Most abundant bacterial genera over time - ", p),
         subtitle = "Top 8 genera by mean Bracken reads (zero-filled) - dominant members, not selected for variability",
         x = "Collection date", y = "Estimated reads (Bracken)") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold.italic"))
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 10.2. VIRAL - trajectories of the most ABUNDANT vOTUs (by mean CoverM
#       coverage), per participant - viral counterpart to 10.1.
# -----------------------------------------------------------------------

top_abundant_votus <- votu_reads_full |>
  group_by(participant, contig) |>
  summarise(mean_abund = mean(abund), .groups = "drop") |>
  group_by(participant) |>
  slice_max(mean_abund, n = 8) |>
  arrange(participant, desc(mean_abund))

message("=== Most abundant vOTUs per participant (by mean CoverM coverage, not CV) ===")
top_abundant_votus |>
  left_join(taxonomy |> select(contig, taxon_label), by = "contig") |>
  print(n = Inf)

top_abundant_votu_line_data <- votu_reads_full |>
  semi_join(top_abundant_votus, by = c("participant", "contig")) |>
  left_join(taxonomy |> select(contig, taxon_label), by = "contig") |>
  left_join(votu_display_lookup |> select(contig, votu_display), by = "contig") |>
  mutate(facet_label = paste0(votu_display, "\n", taxon_label))

for (p in unique(top_abundant_votu_line_data$participant)) {
  
  plot_data <- top_abundant_votu_line_data |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = abund)) +
    geom_line(colour = COL_VIRAL) +
    geom_point(size = 1.8, colour = COL_VIRAL) +
    facet_wrap(~ facet_label, scales = "free_y") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Most abundant vOTUs over time - ", p),
         subtitle = "Top 8 vOTUs by mean CoverM trimmed-mean coverage - dominant members, not selected for variability",
         x = "Collection date", y = "CoverM trimmed-mean coverage") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 10.3. BACTERIAL - single heatmap view: top 20 genera x every sample,
#       per participant - the "what's there and how much" picture in one
#       image, rather than reading it off 8 separate line charts.
# -----------------------------------------------------------------------

top20_genera_by_participant <- genus_reads_full |>
  group_by(participant, bacterial_genus) |>
  summarise(mean_reads = mean(reads), .groups = "drop") |>
  group_by(participant) |>
  slice_max(mean_reads, n = 20) |>
  ungroup()

genus_heatmap_data <- genus_reads_full |>
  semi_join(top20_genera_by_participant, by = c("participant", "bacterial_genus")) |>
  group_by(participant) |>
  mutate(bacterial_genus = fct_reorder(bacterial_genus, reads, .fun = mean)) |>
  ungroup()

for (p in unique(genus_heatmap_data$participant)) {
  
  plot_data <- genus_heatmap_data |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = bacterial_genus, fill = log1p(reads))) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_viridis_c(option = "magma", name = "log(reads + 1)") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Bacterial genus abundance heatmap - ", p),
         subtitle = "Top 20 genera by mean abundance, every sample (zero-filled) - log-scaled for visibility across orders of magnitude",
         x = "Collection date", y = NULL) +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8, face = "italic"),
          legend.position = "right")
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 10.3B. BACTERIAL - single heatmap view: top 20 SPECIES x every sample,
#        per participant (NEW 2026-09-30) - species-level counterpart to
#        10.3 above, same "what's there and how much" picture but at
#        finer taxonomic resolution. Reuses species_reads_full (Part 3,
#        Section 3.6C) - already zero-filled across each participant's
#        full sample set and already carries collection_date, same as
#        genus_reads_full above. Kept as its own section (not merged
#        into 10.3) so a reader can go to genus OR species directly
#        without paging past the other.
# -----------------------------------------------------------------------

top20_species_by_participant <- species_reads_full |>
  group_by(participant, bacterial_taxon) |>
  summarise(mean_reads = mean(reads), .groups = "drop") |>
  group_by(participant) |>
  slice_max(mean_reads, n = 20) |>
  ungroup()

species_heatmap_data <- species_reads_full |>
  semi_join(top20_species_by_participant, by = c("participant", "bacterial_taxon")) |>
  group_by(participant) |>
  mutate(bacterial_taxon = fct_reorder(bacterial_taxon, reads, .fun = mean)) |>
  ungroup()

for (p in unique(species_heatmap_data$participant)) {
  
  plot_data <- species_heatmap_data |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = bacterial_taxon, fill = log1p(reads))) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_viridis_c(option = "magma", name = "log(reads + 1)") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Bacterial species abundance heatmap - ", p),
         subtitle = "Top 20 species by mean abundance, every sample (zero-filled) - log-scaled for visibility across orders of magnitude",
         x = "Collection date", y = NULL) +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7, face = "italic"),
          legend.position = "right")
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}

# -----------------------------------------------------------------------
# 10.4. VIRAL - single heatmap view: top 20 vOTUs x every sample, per
#       participant - viral counterpart to 10.3.
# -----------------------------------------------------------------------

top20_votus_by_participant <- votu_reads_full |>
  group_by(participant, contig) |>
  summarise(mean_abund = mean(abund), .groups = "drop") |>
  group_by(participant) |>
  slice_max(mean_abund, n = 20) |>
  ungroup()

votu_heatmap_data <- votu_reads_full |>
  semi_join(top20_votus_by_participant, by = c("participant", "contig")) |>
  left_join(taxonomy |> select(contig, taxon_label), by = "contig") |>
  left_join(votu_display_lookup |> select(contig, votu_display), by = "contig") |>
  mutate(facet_label = paste0(votu_display, " (", taxon_label, ")")) |>
  group_by(participant) |>
  mutate(facet_label = fct_reorder(facet_label, abund, .fun = mean)) |>
  ungroup()

for (p in unique(votu_heatmap_data$participant)) {
  
  plot_data <- votu_heatmap_data |> filter(participant == p)
  
  gg <- ggplot(plot_data, aes(x = collection_date, y = facet_label, fill = log1p(abund))) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_viridis_c(option = "viridis", name = "log(coverage + 1)") +
    scale_x_date(breaks = sort(unique(plot_data$collection_date)), date_labels = "%d %b %y") +
    labs(title = paste0("Viral vOTU abundance heatmap - ", p),
         subtitle = "Top 20 vOTUs by mean CoverM coverage, every sample - log-scaled for visibility across orders of magnitude",
         x = "Collection date", y = NULL) +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7),
          legend.position = "right")
  
  if (p == "L1") gg <- gg + antibiotic_marker_date()
  
  print(style_figure(gg))
}


# #############################################################################
# PART 11: ADDITIONAL ANALYSES (NEW 2026-10-07)
# Seven add-ons to Aims A-C, all built from objects that already exist by
# this point in the script (viral_long, bacterial_long, abundance_wide,
# genus_reads_full, sample_dates, Part 5's clr_transform_matrix(), Part 8's
# get_participant_clr(), Part 9's window_wilcoxon_test()):
#   11.1  Ordination (PCoA of Aitchison distance)            -> Aim A / C
#   11.2  Within-person vs between-person distance           -> Aim A
#   11.3  Consecutive-sample membership turnover             -> Aim A / C
#   11.4  Sequencing-depth sanity check on richness          -> QC / Aim A
#   11.4B Depth follow-up: L1 depth pre/post, distance vs depth,
#         rarefied distance, depth-adjusted richness           -> QC / Aim A, C
#   11.5  Which species (11.5B) and genera (11.5) responded
#         to L1's antibiotic                                  -> Aim C
#   11.6  Temperate vs virulent fraction pre/post antibiotic -> Aim B / C
#   11.7  Predicted host genus breakdown                     -> Aim C
# Everything here is exploratory/descriptive unless a test is named - at
# 3 participants and 4-5 samples per antibiotic window, the point is
# showing structure honestly, not squeezing out significance.
# #############################################################################

# -----------------------------------------------------------------------
# 11.0. Shared setup: sample metadata + one CLR matrix per domain across
#       ALL participants (taxa x samples, zero-filled, then Part 5's
#       clr_transform_matrix() - same treatment as Parts 5/8/9).
# -----------------------------------------------------------------------

sample_meta <- viral_long |>
  distinct(sample_id, participant, timepoint, collection_date) |>
  mutate(phase = case_when(
    participant != "L1" ~ "no perturbation",
    timepoint <= 4      ~ "L1 pre-antibiotic",
    TRUE                ~ "L1 post-antibiotic"
  )) |>
  arrange(participant, timepoint)

get_combined_clr <- function(domain = c("viral", "bacterial")) {
  domain <- match.arg(domain)
  ids <- sample_meta$sample_id
  
  if (domain == "viral") {
    mat <- abundance_wide |>
      select(contig, all_of(ids)) |>
      column_to_rownames("contig") |>
      as.matrix()
  } else {
    mat <- bacterial_long |>
      filter(sample_id %in% ids) |>
      select(bacterial_taxon, sample_id, new_est_reads) |>
      pivot_wider(names_from = sample_id, values_from = new_est_reads,
                  values_fill = 0, values_fn = sum) |>
      column_to_rownames("bacterial_taxon") |>
      as.matrix()
    mat <- mat[, intersect(ids, colnames(mat)), drop = FALSE]
  }
  
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  clr_transform_matrix(mat)
}

combined_clr <- list(viral = get_combined_clr("viral"),
                     bacterial = get_combined_clr("bacterial"))

# -----------------------------------------------------------------------
# 11.1. ORDINATION - PCoA of Aitchison distance (Euclidean on CLR)
# -----------------------------------------------------------------------
# Why: Part 8 reports L1's post-antibiotic shift as a number (+77%
# Aitchison distance from baseline); an ordination shows where each sample
# actually sits relative to the others. With only 3 people, between-person
# differences will dominate axes 1-2 of the all-participant plot, so a
# second L1-only ordination (11.1b) is what actually shows the antibiotic
# shift. Base-R cmdscale() on the Euclidean CLR distance (= PCA of the CLR
# matrix) - no extra package needed.

run_pcoa <- function(clr_mat) {
  d   <- dist(t(clr_mat))
  ord <- cmdscale(d, k = 2, eig = TRUE)
  pct <- 100 * ord$eig[1:2] / sum(ord$eig[ord$eig > 0])
  list(
    scores = tibble(sample_id = rownames(ord$points),
                    PC1 = ord$points[, 1], PC2 = ord$points[, 2]),
    pct = pct
  )
}

ordination_data <- purrr::map_dfr(c("viral", "bacterial"), function(dom) {
  o <- run_pcoa(combined_clr[[dom]])
  o$scores |>
    left_join(sample_meta, by = "sample_id") |>
    mutate(panel = sprintf("%s (PC1 %.1f%%, PC2 %.1f%%)",
                           str_to_title(dom), o$pct[1], o$pct[2]))
}) |>
  mutate(panel = factor(panel, levels = unique(panel)))

# Print the axis percentages (and make them quotable without opening the figure)
message("=== PCoA % variance explained, all participants (also in the panel titles) ===")
print(ordination_data |> distinct(panel))

p_ordination <- ggplot(ordination_data,
                       aes(x = PC1, y = PC2, colour = participant, shape = phase)) +
  geom_path(aes(group = participant), alpha = 0.35, linewidth = 0.4) +
  geom_point(size = 3) +
  geom_text(aes(label = timepoint), vjust = -1, size = 2.6, show.legend = FALSE) +
  facet_wrap(~ panel, scales = "free") +
  scale_colour_manual(values = PARTICIPANT_COLORS) +
  scale_shape_manual(values = c("no perturbation" = 16,
                                "L1 pre-antibiotic" = 1,
                                "L1 post-antibiotic" = 17)) +
  labs(title = "Sample ordination (PCoA of Aitchison distance), all participants",
       subtitle = "Each point is one sample, numbered by timepoint and joined in time order within participant. Hollow circle = L1 before antibiotic, triangle = L1 after.",
       x = "PC1", y = "PC2", colour = "Participant", shape = NULL) +
  theme_thesis()

print(style_figure(p_ordination, width = 11, height = 5.5))

# 11.1b. L1 only - same ordination, restricted to L1's own samples, so the
# antibiotic contrast isn't swamped by between-person variation.
l1_ordination_data <- purrr::map_dfr(c("viral", "bacterial"), function(dom) {
  o <- run_pcoa(get_participant_clr("L1", dom))
  o$scores |>
    left_join(sample_meta, by = "sample_id") |>
    mutate(panel = sprintf("%s (PC1 %.1f%%, PC2 %.1f%%)",
                           str_to_title(dom), o$pct[1], o$pct[2]))
}) |>
  mutate(panel = factor(panel, levels = unique(panel)))

message("=== PCoA % variance explained, L1 only (also in the panel titles) ===")
print(l1_ordination_data |> distinct(panel))
message("=== L1 ordination: pre- vs post-antibiotic samples along PC1 and PC2 (range overlap = do the groups separate?) ===")
print(l1_ordination_data |>
        group_by(panel, phase) |>
        summarise(n = n(),
                  PC1_min = min(PC1), PC1_max = max(PC1),
                  PC2_min = min(PC2), PC2_max = max(PC2), .groups = "drop") |>
        mutate(across(where(is.double), ~ round(.x, 2))))

p_ordination_l1 <- ggplot(l1_ordination_data,
                          aes(x = PC1, y = PC2, colour = phase)) +
  geom_path(aes(group = 1), colour = "grey60", linewidth = 0.4) +
  geom_point(size = 3.2) +
  geom_text(aes(label = timepoint), vjust = -1, size = 2.8, show.legend = FALSE) +
  facet_wrap(~ panel, scales = "free") +
  scale_colour_manual(values = c("L1 pre-antibiotic" = COL_DIV_HIGH,
                                 "L1 post-antibiotic" = COL_DIV_LOW)) +
  labs(title = "L1 sample ordination (PCoA of Aitchison distance), before vs after antibiotic",
       subtitle = "L1's own samples only, numbered by timepoint and joined in time order. Axes are recomputed within L1, so they are not comparable with the all-participant ordination.",
       x = "PC1", y = "PC2", colour = NULL) +
  theme_thesis()

print(style_figure(p_ordination_l1, width = 11, height = 5.5))

# -----------------------------------------------------------------------
# 11.2. WITHIN-PERSON vs BETWEEN-PERSON DISTANCE
# -----------------------------------------------------------------------
# Why: the standard stability result in longitudinal microbiome work - is a
# person more similar to their own other samples than to other people's?
# DESCRIPTIVE ONLY: pairwise distances share samples, so they aren't
# independent and a p-value across them would be misleading; report the
# medians and the between:within ratio instead. "Within L1" includes
# pre-vs-post antibiotic pairs, so it is inflated relative to D3/G2 by
# design - that is itself part of the L1 story, not a flaw.

pair_distances <- purrr::map_dfr(c("viral", "bacterial"), function(dom) {
  dm  <- as.matrix(dist(t(combined_clr[[dom]])))
  ids <- rownames(dm)
  idx <- which(upper.tri(dm), arr.ind = TRUE)
  tibble(domain = dom, sample_a = ids[idx[, 1]], sample_b = ids[idx[, 2]],
         distance = dm[idx])
}) |>
  left_join(sample_meta |> select(sample_id, participant_a = participant),
            by = c("sample_a" = "sample_id")) |>
  left_join(sample_meta |> select(sample_id, participant_b = participant),
            by = c("sample_b" = "sample_id")) |>
  mutate(
    type = if_else(participant_a == participant_b, "within", "between"),
    comparison = if_else(type == "within", paste0("Within ", participant_a),
                         "Between participants"),
    comparison = factor(comparison,
                        levels = c(paste0("Within ", sort(unique(sample_meta$participant))),
                                   "Between participants"))
  )

within_between_summary <- pair_distances |>
  group_by(domain, comparison) |>
  summarise(n_pairs = n(), median_distance = median(distance),
            mean_distance = mean(distance), .groups = "drop")

within_between_ratio <- pair_distances |>
  group_by(domain, type) |>
  summarise(median_distance = median(distance), .groups = "drop") |>
  pivot_wider(names_from = type, values_from = median_distance) |>
  mutate(ratio_between_to_within = between / within)

message("=== Pairwise Aitchison distance: within-person vs between-person (descriptive; pairs are not independent) ===")
print(within_between_summary, n = Inf)
print(within_between_ratio)

p_within_between <- ggplot(pair_distances, aes(x = comparison, y = distance, fill = comparison)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.3) +
  facet_wrap(~ domain, scales = "free_y", labeller = labeller(domain = str_to_title)) +
  scale_fill_manual(values = c(setNames(PARTICIPANT_COLORS, paste0("Within ", names(PARTICIPANT_COLORS))),
                               "Between participants" = "grey70"),
                    guide = "none") +
  labs(title = "Within-person vs between-person community distance (Aitchison)",
       subtitle = "Every pair of samples; 'Within L1' includes pre-vs-post-antibiotic pairs. Descriptive only - pairwise distances share samples and are not independent.",
       x = NULL, y = "Aitchison distance") +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(style_figure(p_within_between, width = 10, height = 5.5))

# -----------------------------------------------------------------------
# 11.3. CONSECUTIVE-SAMPLE MEMBERSHIP TURNOVER
# -----------------------------------------------------------------------
# Why: Aim A's persistence claim rests on detection frequency + a KS test;
# this measures turnover directly - how much of the detected membership is
# different between each sample and the next one - viral (vOTUs) vs
# bacterial (species), per participant. Jaccard dissimilarity = 1 -
# shared/union of detected features. Sparse samples inflate it for both
# domains, so read viral vs bacterial within a participant, not absolute
# values across them. Gaps between visits are irregular (days_between is
# kept in the table), so it is not rescaled to a per-day rate.

presence_long <- bind_rows(
  viral_long |> filter(abundance > 0) |>
    transmute(domain = "viral", participant, sample_id, feature = contig),
  bacterial_long |> filter(new_est_reads > 0) |>
    transmute(domain = "bacterial", participant, sample_id, feature = bacterial_taxon)
)

consecutive_turnover <- purrr::map_dfr(c("viral", "bacterial"), function(dom) {
  purrr::map_dfr(all_participants, function(part) {
    meta <- sample_meta |> filter(participant == part) |> arrange(timepoint)
    pres <- presence_long |> filter(domain == dom, participant == part)
    sets <- split(pres$feature, pres$sample_id)
    
    purrr::map_dfr(seq_len(nrow(meta) - 1), function(i) {
      a <- sets[[meta$sample_id[i]]]
      b <- sets[[meta$sample_id[i + 1]]]
      if (is.null(a)) a <- character(0)
      if (is.null(b)) b <- character(0)
      union_n <- length(union(a, b))
      tibble(
        domain = dom, participant = part,
        from_sample = meta$sample_id[i], to_sample = meta$sample_id[i + 1],
        to_date = meta$collection_date[i + 1],
        days_between = as.numeric(meta$collection_date[i + 1] - meta$collection_date[i]),
        jaccard_dissimilarity = if (union_n == 0) NA_real_
        else 1 - length(intersect(a, b)) / union_n,
        crosses_antibiotic = part == "L1" && meta$timepoint[i] == 4 && meta$timepoint[i + 1] == 5
      )
    })
  })
})

turnover_domain_test <- consecutive_turnover |>
  select(participant, from_sample, domain, jaccard_dissimilarity) |>
  pivot_wider(names_from = domain, values_from = jaccard_dissimilarity) |>
  group_by(participant) |>
  summarise(
    n_transitions    = n(),
    median_viral     = median(viral, na.rm = TRUE),
    median_bacterial = median(bacterial, na.rm = TRUE),
    wilcoxon_paired_p = suppressWarnings(tryCatch(
      wilcox.test(viral, bacterial, paired = TRUE)$p.value,
      error = function(e) NA_real_)),
    .groups = "drop"
  ) |>
  mutate(p_adj_BH = p.adjust(wilcoxon_paired_p, method = "BH"))

message("=== Consecutive-sample turnover (Jaccard dissimilarity of detected membership): viral vs bacterial, paired Wilcoxon per participant ===")
print(turnover_domain_test)
message("CAVEAT: 7-9 transitions per participant, so the paired test has limited power; the pattern across all three participants matters more than any single p-value.")

p_turnover <- ggplot(consecutive_turnover,
                     aes(x = to_date, y = jaccard_dissimilarity, colour = domain)) +
  geom_line() +
  geom_point(size = 2) +
  antibiotic_marker_date() +
  facet_wrap(~ participant, scales = "free_x") +
  scale_colour_manual(values = c(viral = COL_VIRAL, bacterial = COL_BACTERIAL)) +
  scale_y_continuous(limits = c(0, 1)) +
  participant_date_scales() +
  labs(title = "Turnover between consecutive samples (Jaccard dissimilarity of detected membership)",
       subtitle = "Each point is the change in detected vOTUs (viral) or species (bacterial) between a sample and the one before it, plotted at the later sample's date. 1 = nothing shared, 0 = identical.",
       x = "Collection date (later sample of each pair)", y = "Jaccard dissimilarity", colour = NULL) +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(style_figure(p_turnover, width = 10, height = 5))

# -----------------------------------------------------------------------
# 11.4. SEQUENCING-DEPTH SANITY CHECK
# -----------------------------------------------------------------------
# Why: if richness simply tracks how many reads a sample has, some of the
# diversity fluctuation in Aim A is a sequencing artefact. Depth proxy =
# total Bracken-assigned reads per sample (after phage/human filtering) -
# the only per-sample read total this script has access to - so it
# measures depth of the bacterial classification, used here as a proxy for
# overall library size for both domains. Spearman correlation, since the
# relationship needn't be linear.

sample_depth <- bacterial_long |>
  group_by(sample_id) |>
  summarise(bracken_reads = sum(new_est_reads), .groups = "drop")

depth_richness <- bind_rows(
  viral_long |> group_by(sample_id) |>
    summarise(richness = sum(abundance > 0), .groups = "drop") |>
    mutate(domain = "viral"),
  bacterial_long |> group_by(sample_id) |>
    summarise(richness = n_distinct(bacterial_taxon[new_est_reads > 0]), .groups = "drop") |>
    mutate(domain = "bacterial")
) |>
  left_join(sample_depth, by = "sample_id") |>
  left_join(sample_meta |> select(sample_id, participant, timepoint), by = "sample_id")

depth_cor_row <- function(d) {
  tibble(n_samples = nrow(d),
         spearman_rho = suppressWarnings(cor(d$bracken_reads, d$richness, method = "spearman")),
         spearman_p = suppressWarnings(tryCatch(
           cor.test(d$bracken_reads, d$richness, method = "spearman")$p.value,
           error = function(e) NA_real_)))
}

depth_correlations <- bind_rows(
  depth_richness |> group_by(domain) |> group_modify(~ depth_cor_row(.x)) |>
    ungroup() |> mutate(scope = "all samples"),
  depth_richness |> group_by(domain, participant) |> group_modify(~ depth_cor_row(.x)) |>
    ungroup() |> mutate(scope = participant)
) |>
  select(domain, scope, n_samples, spearman_rho, spearman_p)

message("=== Richness vs sequencing depth (Bracken-assigned reads as depth proxy), Spearman ===")
print(depth_correlations, n = Inf)
message("READING THIS: a strong positive rho means richness partly tracks depth - report diversity results with that caveat. A weak/absent rho supports treating the richness fluctuation as biological.")

p_depth <- ggplot(depth_richness, aes(x = bracken_reads, y = richness, colour = participant)) +
  geom_point(size = 2.4) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey40",
              linewidth = 0.5, linetype = "dashed") +
  facet_wrap(~ domain, scales = "free_y", labeller = labeller(domain = str_to_title)) +
  scale_x_log10(labels = scales::comma) +
  scale_colour_manual(values = PARTICIPANT_COLORS) +
  labs(title = "Richness vs sequencing depth",
       subtitle = "Each point is one sample. Depth proxy = total Bracken-assigned reads (log scale); dashed line = linear fit across all samples, for orientation only.",
       x = "Bracken-assigned reads per sample (depth proxy)", y = "Richness (features detected)",
       colour = "Participant") +
  theme_thesis()

print(style_figure(p_depth, width = 10, height = 5))

# -----------------------------------------------------------------------
# 11.4B. DEPTH FOLLOW-UP (added 2026-10-07)
# -----------------------------------------------------------------------
# 11.4 showed richness tracks depth (rho ~0.8). Three questions follow:
#   (i)   Did L1's depth itself change across the antibiotic? (If post-
#         antibiotic samples are simply deeper/shallower, richness and
#         distance-from-baseline shifts could be depth, not biology.)
#   (ii)  Does distance-from-baseline (Part 8) track depth, and does L1's
#         post-vs-pre difference survive adjusting for it? For bacteria,
#         also a rarefaction robustness check on the Aitchison distance.
#   (iii) Do the Part 9 richness trend/window results change when richness
#         is depth-adjusted?
# Depth proxy is still total Bracken-assigned reads (log10). VIRAL richness
# cannot be rarefied - there is no per-sample viral read total here - so it
# is only adjusted by regression on the bacterial depth proxy, which assumes
# bacterial and viral library sizes move together. Treat viral adjustments
# as approximate. All exploratory: 9 (L1/D3) or fewer samples per person.

# --- (i) depth before vs after the antibiotic ---------------------------
depth_by_sample <- sample_depth |>
  left_join(sample_meta |> select(sample_id, participant, timepoint, collection_date),
            by = "sample_id") |>
  mutate(window = if_else(timepoint <= 4, "pre", "post"),
         log10_depth = log10(bracken_reads))

depth_prepost <- depth_by_sample |>
  group_by(participant) |>
  group_modify(function(d, key) {
    wt <- window_wilcoxon_test(transmute(d, timepoint, value = bracken_reads))
    tibble(median_reads_pre  = median(d$bracken_reads[d$window == "pre"]),
           median_reads_post = median(d$bracken_reads[d$window == "post"]),
           post_over_pre     = median_reads_post / median_reads_pre,
           n_pre = wt$n_pre, n_post = wt$n_post, wilcoxon_p = wt$wilcoxon_p)
  }) |>
  ungroup()

message("=== (i) Sequencing depth (Bracken-assigned reads) pre vs post, per participant (L1 = real antibiotic; D3/G2 = matched reference) ===")
print(depth_prepost, n = Inf)
message("READING THIS: if L1's post/pre depth ratio is near 1 and no more extreme than D3/G2's, depth does not explain L1's pre/post shifts; a large ratio would mean it could.")

p_depth_time <- ggplot(depth_by_sample,
                       aes(x = collection_date, y = bracken_reads, colour = participant)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.2) +
  antibiotic_marker_date() +
  facet_wrap(~ participant, scales = "free_x") +
  scale_y_log10(labels = scales::comma) +
  scale_colour_manual(values = PARTICIPANT_COLORS, guide = "none") +
  participant_date_scales() +
  labs(title = "Sequencing depth over time",
       subtitle = "Total Bracken-assigned reads per sample (log scale); dashed line = L1's antibiotic course. D3/G2 have no event - shown as reference.",
       x = "Collection date", y = "Bracken-assigned reads (depth proxy)") +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(style_figure(p_depth_time, width = 10, height = 4.5))

# --- (ii) distance-from-baseline vs depth --------------------------------
dist_depth <- resilience_trajectory |>
  left_join(depth_by_sample |> select(sample_id, log10_depth, window), by = "sample_id")

dist_depth_summary <- dist_depth |>
  group_by(domain, participant) |>
  group_modify(function(d, key) {
    rho <- suppressWarnings(cor(d$log10_depth, d$aitchison_distance, method = "spearman"))
    rho_p <- suppressWarnings(tryCatch(
      cor.test(d$log10_depth, d$aitchison_distance, method = "spearman")$p.value,
      error = function(e) NA_real_))
    d$window <- factor(d$window, levels = c("pre", "post"))
    raw_fit <- lm(aitchison_distance ~ window, data = d)
    adj_fit <- lm(aitchison_distance ~ window + log10_depth, data = d)
    cf_raw <- summary(raw_fit)$coefficients
    cf_adj <- summary(adj_fit)$coefficients
    tibble(n = nrow(d),
           spearman_rho_dist_vs_depth = rho, spearman_p = rho_p,
           post_minus_pre_raw      = unname(cf_raw["windowpost", "Estimate"]),
           post_minus_pre_raw_p    = unname(cf_raw["windowpost", "Pr(>|t|)"]),
           post_minus_pre_adj      = unname(cf_adj["windowpost", "Estimate"]),
           post_minus_pre_adj_p    = unname(cf_adj["windowpost", "Pr(>|t|)"]))
  }) |>
  ungroup()

message("=== (ii) Aitchison distance-from-baseline vs depth, and post-minus-pre difference before/after adjusting for log10 depth ===")
print(dist_depth_summary, n = Inf)
message("READING THIS: 'post_minus_pre_raw' is the plain pre/post difference in distance; '_adj' is the same after adding log10 depth to the model. If L1's adjusted difference stays close to its raw value (and larger than D3/G2's), the shift is not a depth artefact. Note baseline samples sit nearest the baseline centroid by construction, and the p-values come from 9 or fewer points - descriptive only.")

# Bacterial rarefaction robustness: resample every sample to the same depth
# (multinomial draw from its own composition), recompute Aitchison distance
# from the pre-antibiotic centroid, average over replicates.
set.seed(20261007)
N_RAREFY_REPS <- 20

bacterial_counts_by_participant <- function(part) {
  bacterial_long |>
    filter(participant == part) |>
    select(bacterial_taxon, sample_id, new_est_reads) |>
    pivot_wider(names_from = sample_id, values_from = new_est_reads, values_fill = 0) |>
    column_to_rownames("bacterial_taxon") |>
    as.matrix()
}

rarefy_depth <- floor(min(sample_depth$bracken_reads))

rarefied_distance_one <- function(part) {
  mat <- bacterial_counts_by_participant(part)
  pre_ids <- sample_meta$sample_id[sample_meta$participant == part & sample_meta$timepoint <= 4]
  reps <- lapply(seq_len(N_RAREFY_REPS), function(r) {
    rar <- apply(mat, 2, function(x) as.numeric(rmultinom(1, size = rarefy_depth, prob = x / sum(x))))
    rownames(rar) <- rownames(mat)
    rar <- rar[rowSums(rar) > 0, , drop = FALSE]
    distance_from_baseline(clr_transform_matrix(rar), pre_ids)
  })
  bind_rows(reps) |>
    group_by(sample_id) |>
    summarise(rarefied_distance = mean(aitchison_distance), .groups = "drop") |>
    mutate(participant = part)
}

rarefied_distances <- purrr::map_dfr(all_participants, rarefied_distance_one) |>
  left_join(sample_meta |> select(sample_id, timepoint), by = "sample_id") |>
  left_join(resilience_trajectory |> filter(domain == "bacterial") |>
              select(sample_id, unrarefied_distance = aitchison_distance), by = "sample_id") |>
  mutate(window = if_else(timepoint <= 4, "pre", "post"))

rarefied_summary <- rarefied_distances |>
  group_by(participant, window) |>
  summarise(unrarefied = median(unrarefied_distance),
            rarefied   = median(rarefied_distance), .groups = "drop") |>
  pivot_wider(names_from = window, values_from = c(unrarefied, rarefied)) |>
  mutate(pct_change_unrarefied = 100 * (unrarefied_post / unrarefied_pre - 1),
         pct_change_rarefied   = 100 * (rarefied_post / rarefied_pre - 1))

message(sprintf("=== (ii-b) Bacterial Aitchison distance-from-baseline, unrarefied vs rarefied to %s reads (%d multinomial replicates, median per window) ===",
                format(rarefy_depth, big.mark = ","), N_RAREFY_REPS))
print(rarefied_summary, n = Inf)
message("READING THIS: if L1's pct_change_rarefied stays close to pct_change_unrarefied (earlier: +77%), the bacterial disruption survives equalising depth. Rarefaction to the shallowest sample discards information and shrinks absolute distances, so compare the percent change, not the raw distances.")

# --- (iii) depth-adjusted richness ---------------------------------------
# Residual richness: richness minus its expectation from log10 depth (one
# regression per domain across all samples). Plus, for bacteria, Hurlbert's
# expected richness at the shallowest sample's depth (rarefied richness).

hurlbert_richness <- function(counts, n) {
  counts <- round(counts[counts > 0])
  N <- sum(counts)
  sum(1 - exp(lchoose(N - counts, n) - lchoose(N, n)))
}

bacterial_rarefied_richness <- bacterial_long |>
  group_by(sample_id) |>
  summarise(value = hurlbert_richness(new_est_reads, rarefy_depth), .groups = "drop") |>
  left_join(sample_meta |> select(sample_id, participant, timepoint), by = "sample_id") |>
  mutate(domain = "bacterial", version = "rarefied (Hurlbert)")

richness_adjusted <- depth_richness |>
  group_by(domain) |>
  mutate(value = resid(lm(richness ~ log10(bracken_reads))),
         version = "depth-residual") |>
  ungroup() |>
  select(sample_id, participant, timepoint, domain, value, version)

richness_raw <- depth_richness |>
  transmute(sample_id, participant, timepoint, domain, value = richness, version = "raw")

richness_versions <- bind_rows(richness_raw, richness_adjusted,
                               bacterial_rarefied_richness |>
                                 select(sample_id, participant, timepoint, domain, value, version))

richness_depth_tests <- richness_versions |>
  group_by(domain, version, participant) |>
  group_modify(function(d, key) {
    kt <- suppressWarnings(tryCatch(cor.test(d$timepoint, d$value, method = "kendall"),
                                    error = function(e) NULL))
    wt <- window_wilcoxon_test(d)
    tibble(kendall_tau = if (is.null(kt)) NA_real_ else unname(kt$estimate),
           kendall_p   = if (is.null(kt)) NA_real_ else kt$p.value,
           n_pre = wt$n_pre, n_post = wt$n_post, wilcoxon_p = wt$wilcoxon_p)
  }) |>
  ungroup() |>
  group_by(version) |>
  # BH pooled across the participants and both domains within each version
  # (six tests for raw / depth-residual, three for rarefied). The raw rows'
  # headline adjusted p-values are the Part 2.6 ones (12-test family,
  # richness + Shannon), joined in below as kendall_p_adj_BH_part26.
  mutate(kendall_p_adj_BH_within_version = p.adjust(kendall_p, method = "BH")) |>
  ungroup() |>
  mutate(version = factor(version, levels = c("raw", "depth-residual", "rarefied (Hurlbert)"))) |>
  arrange(domain, participant, version)

# Cross-check: raw rows here must reproduce Part 9's richness results.
raw_check <- richness_depth_tests |>
  filter(version == "raw") |>
  mutate(domain_cap = str_to_title(domain)) |>
  left_join(diversity_trend_tests |> filter(metric == "richness") |>
              select(participant, domain_cap = domain, part9_kendall_p = kendall_p),
            by = c("participant", "domain_cap"))
message(sprintf("Cross-check vs Part 9 raw richness Kendall p-values: max abs difference = %.2e (should be ~0)",
                max(abs(raw_check$kendall_p - raw_check$part9_kendall_p), na.rm = TRUE)))

richness_depth_tests <- richness_depth_tests |>
  left_join(diversity_trend_tests |>
              filter(metric == "richness") |>
              transmute(participant, domain = str_to_lower(domain), version = "raw",
                        kendall_p_adj_BH_part26 = kendall_p_adj_BH) |>
              mutate(version = factor(version, levels = levels(richness_depth_tests$version))),
            by = c("participant", "domain", "version"))

message("=== (iii) Richness trend (Kendall vs timepoint) and pre/post Wilcoxon: raw vs depth-adjusted ===")
message("kendall_p_adj_BH_part26 = the Part 2.6 12-test BH value (raw rows only; this is the one to report). kendall_p_adj_BH_within_version = BH across participants and domains within that version, for comparing the adjusted versions among themselves.")
print(richness_depth_tests |> select(domain, participant, version, kendall_tau, kendall_p,
                                     kendall_p_adj_BH_part26, kendall_p_adj_BH_within_version,
                                     n_pre, n_post, wilcoxon_p), n = Inf)
message("READING THIS: compare each participant's raw row with its depth-adjusted rows. A trend or pre/post difference that disappears after adjustment was depth; one that persists is not explained by it. Viral rows use the bacterial depth proxy and cannot be rarefied. Same low-power caveat as Part 9 (4 vs 5 samples; smallest attainable Wilcoxon p = 0.016).")

# -----------------------------------------------------------------------
# Detection-count helpers for 11.5 / 11.5B (added 2026-10-07)
# -----------------------------------------------------------------------
# A CLR change for a taxon that is ABSENT in one window is partly set by the
# zero-replacement pseudocount (min nonzero / 2), not by a measured
# abundance. These helpers add, per taxon and participant, how many pre
# (timepoint <=4) and post samples it was detected in, so a "decrease" can
# be read as "present in 4/4 pre samples, 0/5 post" rather than as a number
# on a floor-limited scale.

detection_counts <- function(reads_df, taxon_col, part) {
  d <- reads_df |>
    filter(participant == part) |>
    mutate(window = if_else(timepoint <= 4, "pre", "post"))
  n_pre_s  <- n_distinct(d$sample_id[d$window == "pre"])
  n_post_s <- n_distinct(d$sample_id[d$window == "post"])
  out <- d |>
    group_by(.data[[taxon_col]]) |>
    summarise(
      n_detected       = sum(reads > 0),
      n_detected_pre   = sum(reads > 0 & window == "pre"),
      n_detected_post  = sum(reads > 0 & window == "post"),
      mean_reads_pre   = round(sum(reads[window == "pre"])  / n_pre_s, 1),
      mean_reads_post  = round(sum(reads[window == "post"]) / n_post_s, 1),
      .groups = "drop") |>
    mutate(n_samples_pre = n_pre_s, n_samples_post = n_post_s,
           detected_pre_to_post = sprintf("%d/%d -> %d/%d",
                                          n_detected_pre, n_pre_s,
                                          n_detected_post, n_post_s))
  out
}

response_print_cols <- function(df, taxon_col) {
  df |> select(all_of(c(taxon_col, "delta_clr", "detected_pre_to_post",
                        "mean_reads_pre", "mean_reads_post", "wilcoxon_p_adj_BH")))
}

vanish_emerge_summary <- function(resp, taxon_col, level) {
  resp |>
    group_by(participant) |>
    summarise(
      taxa_tested = n(),
      vanished = sum(n_detected_pre >= 3 & n_detected_post == 0),
      emerged  = sum(n_detected_post >= 3 & n_detected_pre == 0),
      .groups = "drop") |>
    mutate(level = level, .before = 1)
}

# -----------------------------------------------------------------------
# 11.5. WHICH BACTERIAL GENERA RESPONDED TO L1's ANTIBIOTIC?
# -----------------------------------------------------------------------
# Why: Part 8 shows L1's bacterial community shifted (+77% Aitchison
# distance) but not WHICH genera drove it. Per genus: change in mean CLR
# abundance, post (timepoints >4) minus pre (timepoints <=4); D3 and G2
# get the same split as a structurally matched reference (no perturbation),
# so an L1 change can be judged against ordinary pre/post drift. EXPLORATORY:
# with 4 pre and 5 post samples the smallest attainable two-sided Wilcoxon
# p-value is 2/choose(9,4) = 0.016, so after BH correction across ~100
# genera nothing can reach significance - the ranking is descriptive, and
# the p-values are shown only so that limitation is visible, not hidden.
# Genera need >= 3 detections in a participant's samples to be included.

genus_clr_by_participant <- function(part) {
  mat <- genus_reads_full |>
    filter(participant == part) |>
    select(bacterial_genus, sample_id, reads) |>
    pivot_wider(names_from = sample_id, values_from = reads, values_fill = 0) |>
    column_to_rownames("bacterial_genus") |>
    as.matrix()
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  clr_transform_matrix(mat)
}

genus_response <- purrr::map_dfr(all_participants, function(part) {
  clr <- genus_clr_by_participant(part)
  tps <- genus_reads_full |> filter(participant == part) |> distinct(sample_id, timepoint)
  pre  <- tps$sample_id[tps$timepoint <= 4]
  post <- tps$sample_id[tps$timepoint > 4]
  n_det <- detection_counts(genus_reads_full, "bacterial_genus", part)
  
  tibble(
    participant   = part,
    bacterial_genus = rownames(clr),
    mean_clr_pre  = rowMeans(clr[, pre, drop = FALSE]),
    mean_clr_post = rowMeans(clr[, post, drop = FALSE]),
    delta_clr     = mean_clr_post - mean_clr_pre,
    wilcoxon_p    = purrr::map_dbl(seq_len(nrow(clr)), function(i) {
      suppressWarnings(tryCatch(wilcox.test(clr[i, post], clr[i, pre])$p.value,
                                error = function(e) NA_real_))
    })
  ) |>
    left_join(n_det, by = "bacterial_genus") |>
    filter(n_detected >= 3)
}) |>
  group_by(participant) |>
  mutate(wilcoxon_p_adj_BH = p.adjust(wilcoxon_p, method = "BH")) |>
  ungroup()

l1_genus_ranked <- genus_response |>
  filter(participant == "L1") |>
  arrange(desc(delta_clr))

message("=== L1: genera with the largest increase in CLR abundance post-antibiotic (top 10) ===")
print(response_print_cols(l1_genus_ranked |> slice_head(n = 10), "bacterial_genus"), n = Inf)
message("=== L1: genera with the largest decrease in CLR abundance post-antibiotic (bottom 10) ===")
print(response_print_cols(l1_genus_ranked |> slice_tail(n = 10) |> arrange(delta_clr), "bacterial_genus"), n = Inf)
message(sprintf("L1 genera tested: %d; smallest BH-adjusted p = %.3f (smallest attainable raw p with 4 vs 5 samples = 0.016). Descriptive ranking only.",
                nrow(l1_genus_ranked), min(l1_genus_ranked$wilcoxon_p_adj_BH, na.rm = TRUE)))

l1_genus_selected <- bind_rows(slice_head(l1_genus_ranked, n = 10),
                               slice_tail(l1_genus_ranked, n = 10)) |>
  distinct(bacterial_genus, .keep_all = TRUE)
genus_levels <- l1_genus_selected |> arrange(delta_clr) |> pull(bacterial_genus)

l1_genus_plot_data <- l1_genus_selected |>
  mutate(bacterial_genus = factor(bacterial_genus, levels = genus_levels),
         direction = delta_clr > 0)
genus_reference_data <- genus_response |>
  filter(participant != "L1", bacterial_genus %in% genus_levels) |>
  mutate(bacterial_genus = factor(bacterial_genus, levels = genus_levels))

p_genus_response <- ggplot(l1_genus_plot_data, aes(x = delta_clr, y = bacterial_genus)) +
  geom_col(aes(fill = direction), width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_point(data = genus_reference_data, aes(colour = participant),
             shape = 18, size = 3) +
  scale_fill_manual(values = c("TRUE" = COL_DIV_HIGH, "FALSE" = COL_DIV_LOW),
                    labels = c("TRUE" = "L1: higher after antibiotic",
                               "FALSE" = "L1: lower after antibiotic"),
                    name = NULL) +
  scale_colour_manual(values = PARTICIPANT_COLORS[c("D3", "G2")],
                      name = "Reference, same pre/post split (no antibiotic)") +
  labs(title = "L1: bacterial genera with the largest change in abundance after the antibiotic",
       subtitle = "Change in mean CLR abundance, post minus pre, top 10 increases and decreases (genera detected in >=3 L1 samples). Diamonds = the same genus in D3/G2 over the same sample split. Exploratory - 4 vs 5 samples.",
       x = "Change in mean CLR abundance (post - pre)", y = NULL) +
  theme_thesis() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.box = "vertical")

print(style_figure(p_genus_response, width = 9, height = 7.5))

# -----------------------------------------------------------------------
# 11.5B. SAME QUESTION AT SPECIES LEVEL (added 2026-10-07)
# -----------------------------------------------------------------------
# Species is the more informative resolution for an antibiotic response -
# susceptibility is often species-specific, and a genus-level change can
# average a responder together with non-responders. Same method and same
# exploratory caveat as 11.5 (change in mean CLR abundance, post minus
# pre; D3/G2 as matched reference; >=3 detections required). Costs: species
# data is sparser (more zeros, fewer taxa clear the >=3-detection floor)
# and Bracken's species-level calls are less certain than genus-level, so
# read this alongside 11.5, not instead of it. Reuses species_reads_full
# (Part 3, Section 3.6C - already zero-filled).

species_clr_by_participant <- function(part) {
  mat <- species_reads_full |>
    filter(participant == part) |>
    select(bacterial_taxon, sample_id, reads) |>
    pivot_wider(names_from = sample_id, values_from = reads, values_fill = 0) |>
    column_to_rownames("bacterial_taxon") |>
    as.matrix()
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  clr_transform_matrix(mat)
}

species_response <- purrr::map_dfr(all_participants, function(part) {
  clr <- species_clr_by_participant(part)
  tps <- species_reads_full |> filter(participant == part) |> distinct(sample_id, timepoint)
  pre  <- tps$sample_id[tps$timepoint <= 4]
  post <- tps$sample_id[tps$timepoint > 4]
  n_det <- detection_counts(species_reads_full, "bacterial_taxon", part)
  
  tibble(
    participant     = part,
    bacterial_taxon = rownames(clr),
    mean_clr_pre    = rowMeans(clr[, pre, drop = FALSE]),
    mean_clr_post   = rowMeans(clr[, post, drop = FALSE]),
    delta_clr       = mean_clr_post - mean_clr_pre,
    wilcoxon_p      = purrr::map_dbl(seq_len(nrow(clr)), function(i) {
      suppressWarnings(tryCatch(wilcox.test(clr[i, post], clr[i, pre])$p.value,
                                error = function(e) NA_real_))
    })
  ) |>
    left_join(n_det, by = "bacterial_taxon") |>
    filter(n_detected >= 3)
}) |>
  group_by(participant) |>
  mutate(wilcoxon_p_adj_BH = p.adjust(wilcoxon_p, method = "BH")) |>
  ungroup()

l1_species_ranked <- species_response |>
  filter(participant == "L1") |>
  arrange(desc(delta_clr))

message("=== L1: SPECIES with the largest increase in CLR abundance post-antibiotic (top 10) ===")
print(response_print_cols(l1_species_ranked |> slice_head(n = 10), "bacterial_taxon"), n = Inf)
message("=== L1: SPECIES with the largest decrease in CLR abundance post-antibiotic (bottom 10) ===")
print(response_print_cols(l1_species_ranked |> slice_tail(n = 10) |> arrange(delta_clr), "bacterial_taxon"), n = Inf)
message(sprintf("L1 species tested: %d (detected in >=3 L1 samples); smallest BH-adjusted p = %.3f (smallest attainable raw p with 4 vs 5 samples = 0.016). Descriptive ranking only.",
                nrow(l1_species_ranked), min(l1_species_ranked$wilcoxon_p_adj_BH, na.rm = TRUE)))

message("READING THE DETECTION COLUMN: 'detected_pre_to_post' = samples with the taxon detected / samples, pre -> post (e.g. '4/4 -> 0/5' = present in every pre sample, never seen after). delta_clr for such taxa is partly set by the zero-replacement floor, so rely on the detection counts for 'vanished' claims.")

l1_species_vanished <- species_response |>
  filter(participant == "L1", n_detected_pre >= 3, n_detected_post == 0) |>
  arrange(desc(mean_reads_pre))
l1_species_emerged <- species_response |>
  filter(participant == "L1", n_detected_post >= 3, n_detected_pre == 0) |>
  arrange(desc(mean_reads_post))

message("=== L1: species detected in >=3 pre-antibiotic samples and in NONE post (vanished), by pre-antibiotic mean reads ===")
print(response_print_cols(l1_species_vanished, "bacterial_taxon"), n = 30)
message("=== L1: species detected in >=3 post-antibiotic samples and in NONE pre (emerged), by post-antibiotic mean reads ===")
print(response_print_cols(l1_species_emerged, "bacterial_taxon"), n = 30)

vanish_emerge_comparison <- bind_rows(
  vanish_emerge_summary(species_response, "bacterial_taxon", "species"),
  vanish_emerge_summary(genus_response, "bacterial_genus", "genus"))
message("=== Vanished / emerged taxa per participant (L1 = antibiotic; D3/G2 = same pre/post split, no event) ===")
print(vanish_emerge_comparison, n = Inf)
message("READING THIS: L1 having clearly more vanished taxa than D3/G2 supports a real loss of species; similar counts would mean this is ordinary turnover at 4-vs-5 samples. Window sizes differ slightly by participant, so compare counts as rough orders of magnitude.")

l1_species_selected <- bind_rows(slice_head(l1_species_ranked, n = 10),
                                 slice_tail(l1_species_ranked, n = 10)) |>
  distinct(bacterial_taxon, .keep_all = TRUE)
species_levels <- l1_species_selected |> arrange(delta_clr) |> pull(bacterial_taxon)

l1_species_plot_data <- l1_species_selected |>
  mutate(bacterial_taxon = factor(bacterial_taxon, levels = species_levels),
         direction = delta_clr > 0)
species_reference_data <- species_response |>
  filter(participant != "L1", bacterial_taxon %in% species_levels) |>
  mutate(bacterial_taxon = factor(bacterial_taxon, levels = species_levels))

p_species_response <- ggplot(l1_species_plot_data, aes(x = delta_clr, y = bacterial_taxon)) +
  geom_col(aes(fill = direction), width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_point(data = species_reference_data, aes(colour = participant),
             shape = 18, size = 3) +
  scale_fill_manual(values = c("TRUE" = COL_DIV_HIGH, "FALSE" = COL_DIV_LOW),
                    labels = c("TRUE" = "L1: higher after antibiotic",
                               "FALSE" = "L1: lower after antibiotic"),
                    name = NULL) +
  scale_colour_manual(values = PARTICIPANT_COLORS[c("D3", "G2")],
                      name = "Reference, same pre/post split (no antibiotic)") +
  labs(title = "L1: bacterial species with the largest change in abundance after the antibiotic",
       subtitle = "Change in mean CLR abundance, post minus pre, top 10 increases and decreases (species detected in >=3 L1 samples). Diamonds = the same species in D3/G2 over the same sample split; absent where not detected in >=3 of their samples. Exploratory - 4 vs 5 samples.",
       x = "Change in mean CLR abundance (post - pre)", y = NULL) +
  theme_thesis() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.box = "vertical")

print(style_figure(p_species_response, width = 9, height = 7.5))

# -----------------------------------------------------------------------
# 11.6. TEMPERATE vs VIRULENT FRACTION, PRE vs POST ANTIBIOTIC
# -----------------------------------------------------------------------
# Why: antibiotics can induce prophages, so a rise in the temperate share
# of the detected virome after L1's course would be a notable Aim B/C
# result. Uses PhaBOX2 lifestyle calls (TYPE), so it only covers vOTUs
# with a call (printed below) - Part 4's lifestyle comparison relies on the same calls, so it shares
# this coverage limit. Two measures per sample: share of called vOTUs that are
# temperate (presence-based), and share of GMPR-normalized abundance that
# is temperate. Same pre/post Wilcoxon + D3/G2 reference framing as
# Section 9.7 (window_wilcoxon_test() is reused as-is).

message("=== PhaBOX2 lifestyle calls (TYPE) across the vOTU catalogue ===")
print(viral_long |> distinct(contig, TYPE) |> count(TYPE, sort = TRUE))

lifestyle_by_sample <- viral_long |>
  filter(abundance > 0, !is.na(TYPE), TYPE != "-") |>
  mutate(lifestyle_class = case_when(
    str_detect(tolower(TYPE), "temperate") ~ "temperate",
    str_detect(tolower(TYPE), "virulent")  ~ "virulent",
    TRUE ~ NA_character_
  )) |>
  filter(!is.na(lifestyle_class)) |>
  mutate(norm_abund = abundance / size_factor) |>
  group_by(participant, sample_id, timepoint, collection_date) |>
  summarise(
    n_called = n_distinct(contig),
    n_temperate = n_distinct(contig[lifestyle_class == "temperate"]),
    temperate_fraction_votus = n_temperate / n_called,
    temperate_fraction_abundance = sum(norm_abund[lifestyle_class == "temperate"]) / sum(norm_abund),
    .groups = "drop"
  )

message(sprintf("Samples with at least one lifestyle-called vOTU detected: %d of %d; median lifestyle-called vOTUs per sample: %.0f",
                nrow(lifestyle_by_sample), nrow(sample_meta), median(lifestyle_by_sample$n_called)))

temperate_window_tests <- lifestyle_by_sample |>
  select(participant, timepoint, temperate_fraction_votus, temperate_fraction_abundance) |>
  pivot_longer(starts_with("temperate_fraction"), names_to = "metric", values_to = "value") |>
  filter(!is.na(value)) |>
  group_by(participant, metric) |>
  group_modify(~ window_wilcoxon_test(.x)) |>
  ungroup()

message("=== Temperate fraction, pre- vs post-antibiotic-window Wilcoxon (L1 = real event; D3/G2 = reference) ===")
print(temperate_window_tests, n = Inf)
message("CAVEAT: few lifestyle-called vOTUs per sample and 4-5 samples per window - treat as exploratory, read L1 against the D3/G2 rows.")

p_temperate <- ggplot(lifestyle_by_sample,
                      aes(x = collection_date, y = temperate_fraction_abundance)) +
  geom_line(colour = COL_VIRAL) +
  geom_point(size = 2, colour = COL_VIRAL) +
  antibiotic_marker_date() +
  facet_wrap(~ participant, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  participant_date_scales() +
  labs(title = "Temperate share of the detected virome over time",
       subtitle = "GMPR-normalized abundance of temperate vOTUs as a share of all lifestyle-called vOTUs detected in each sample (PhaBOX2 calls only; uncalled vOTUs excluded).",
       x = "Collection date", y = "Temperate share of abundance") +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(style_figure(p_temperate, width = 10, height = 5))

# -----------------------------------------------------------------------
# 11.7. PREDICTED HOST GENUS BREAKDOWN
# -----------------------------------------------------------------------
# Why: context for Part 7's 2-pair phage:host finding - what hosts did
# PhaBOX2 predict for the vOTU catalogue as a whole, and how much of the
# catalogue has a species-level prediction at all? Species-level host
# calls are collapsed to genus for display. These are PREDICTED hosts, not
# confirmed infections (see Part 7 header).

host_breakdown_votu <- viral_long |>
  distinct(contig, Host) |>
  mutate(
    host_class = case_when(
      is.na(Host) | Host == "-"  ~ "No host predicted",
      Host == "species:bacterium" ~ "Generic 'bacterium' only",
      TRUE                        ~ "Species-level host"
    ),
    host_species = if_else(host_class == "Species-level host",
                           str_remove(Host, "^(?:(?:species|genus):)+"), NA_character_),
    # PhaBOX2 sometimes returns an unresolved label such as
    # "Akkermansiaceae bacterium" or "Lachnospiraceae bacterium 5_1_63FAA":
    # the first word is then a FAMILY/ORDER, not a genus. Keep those out of
    # the genus table (they are counted separately below).
    host_higher_rank_only = !is.na(host_species) & word(host_species, 2) %in% "bacterium",
    host_higher_label = if_else(host_higher_rank_only, word(host_species, 1), NA_character_),
    host_genus = if_else(host_higher_rank_only, NA_character_, word(host_species, 1))
  )

message("=== Host-prediction coverage across the vOTU catalogue ===")
print(host_breakdown_votu |> count(host_class, sort = TRUE))

strip_brackets <- function(x) str_remove_all(x, "[\\[\\]]")
bracken_genera_clean <- unique(strip_brackets(bacterial_long$bacterial_genus))

host_genus_table <- host_breakdown_votu |>
  filter(!is.na(host_genus)) |>
  count(host_genus, name = "n_votus_in_catalogue") |>
  mutate(
    # bracketed names ("[Clostridium]") are compared without the brackets
    genus_detected_in_bracken = strip_brackets(host_genus) %in% bracken_genera_clean,
    note = if_else(str_detect(host_genus, "^\\["),
                   "bracketed (polyphyletic) genus name - NCBI may have renamed it in Bracken's database, so a FALSE here is not strong evidence of absence",
                   NA_character_)
  ) |>
  arrange(desc(n_votus_in_catalogue))

message("=== Predicted host genera: vOTUs in catalogue, and whether that genus appears in Bracken at all ===")
print(host_genus_table, n = Inf)

host_higher_rank_table <- host_breakdown_votu |>
  filter(host_higher_rank_only) |>
  count(host_higher_label, name = "n_votus_in_catalogue") |>
  arrange(desc(n_votus_in_catalogue))
message(sprintf("=== Host predictions resolved only to family/order (%d vOTUs; excluded from the genus table above and from the host-genus figure) ===",
                sum(host_higher_rank_table$n_votus_in_catalogue)))
print(host_higher_rank_table, n = Inf)

host_genus_detected <- viral_long |>
  filter(abundance > 0) |>
  distinct(contig, participant) |>
  inner_join(host_breakdown_votu |> filter(!is.na(host_genus)) |> select(contig, host_genus),
             by = "contig") |>
  count(host_genus, participant, name = "n_votus")

if (nrow(host_genus_detected) > 0) {
  host_levels <- host_genus_detected |>
    group_by(host_genus) |>
    summarise(total = sum(n_votus), .groups = "drop") |>
    arrange(total) |>
    pull(host_genus)
  
  p_host_genus <- ggplot(host_genus_detected |>
                           mutate(host_genus = factor(host_genus, levels = host_levels)),
                         aes(x = n_votus, y = host_genus, fill = participant)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = PARTICIPANT_COLORS) +
    scale_x_continuous(breaks = integer_breaks, expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Predicted host genera of the detected vOTUs, by participant",
         subtitle = sprintf("PhaBOX2 species-level host predictions collapsed to genus, counting vOTUs detected in each participant. %d of %d catalogue vOTUs have a species-level prediction; the rest have none or only a generic 'bacterium' call. Predicted, not confirmed, hosts.",
                            sum(host_breakdown_votu$host_class == "Species-level host"),
                            nrow(host_breakdown_votu)),
         x = "vOTUs detected", y = NULL, fill = "Participant") +
    theme_thesis() +
    theme(axis.text.y = element_text(face = "italic"))
  
  print(style_figure(p_host_genus, width = 8, height = 4.5))
} else {
  message("No detected vOTUs have a species-level host prediction - skipping the host-genus figure.")
}


# #############################################################################
# END-OF-SCRIPT: figure legend text + vOTU display-ID lookup (NEW 2026-09-29)
# #############################################################################
# Every figure's original title/subtitle (moved off the image itself by
# style_figure(), see Part 2) is collected here in the order it was drawn,
# so the actual wording is on hand when writing the thesis's real figure
# legends - copy/adapt from this list rather than reconstructing captions
# from the code. Printed as one block at the very end of the script so
# it's the last thing visible after a full source() run.

message(sprintf("=== Figure titles/subtitles captured for thesis figure legends (%d figures) ===",
                length(figure_captions)))
for (i in seq_along(figure_captions)) {
  ttl <- figure_captions[[i]]$title
  sub <- figure_captions[[i]]$subtitle
  message(sprintf("[%02d] %s%s", i,
                  if (is.null(ttl)) "(no title)" else ttl,
                  if (is.null(sub)) "" else paste0(" -- ", gsub("\n", " ", sub))))
}

# The real contig ID behind every vOTU_### label used on every figure
# above (Part 1, Section 1.1f) - worth including as a supplementary
# table in the thesis so a reader can trace a specific vOTU_### back to
# its actual contig/assembly, and worth double-checking here that the
# count matches the 151-vOTU catalogue size.
message(sprintf("=== vOTU display-ID lookup (%d vOTUs) - see votu_display_lookup in the R environment for the full table ===",
                nrow(votu_display_lookup)))
print(votu_display_lookup, n = Inf)
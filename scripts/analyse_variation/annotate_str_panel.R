library(tidyverse)

source("scripts/analyse_variation/str_gt_utils.R")

df_str_panel <- load_str_panel("/Users/maxverbiest/PhD/data/str_panels/tral_and_perf_panel.tsv")

df_segdup_panel <- load_str_panel("/Users/maxverbiest/PhD/data/str_panels/tral_and_perf_panel_segdups.tsv")   # Segmental duplications
df_str_panel <- add_panel_info(df_str_panel, df_segdup_panel, "in_segdup")

df_exon_panel <- load_str_panel("/Users/maxverbiest/PhD/data/str_panels/tral_and_perf_panel_exons.tsv")       # Exonic regions
df_cds_panel <- load_str_panel("/Users/maxverbiest/PhD/data/str_panels/tral_and_perf_panel_CDS.tsv")          # CDS regions

# Add genomic region type, drop exon and CDS columns
df_str_panel <- add_panel_info(df_str_panel, df_exon_panel, "in_exon")
df_str_panel <- add_panel_info(df_str_panel, df_cds_panel, "in_cds")
df_str_panel <- df_str_panel %>% 
  mutate(
    region_type = case_when(
      !in_exon ~ "intron/intergenic",
      in_exon & !in_cds ~ "UTR",
      in_cds ~ "CDS"
    )
  ) %>% 
  select(-c(in_exon, in_cds))


# load tral_perf panel for STRs that have STRs in their flanking regions (50bp)
df_neighbour_panel <- load_str_panel("/Users/maxverbiest/PhD/data/str_panels/tral_and_perf_panel_flanks.tsv")            

# load tral_perf panel for STRs that have STRs in their flanking regions with exactly the same consensus unit
df_neighbour_matching_panel <- load_str_panel("/Users/maxverbiest/PhD/data/str_panels/tral_and_perf_panel_flanks_cons.tsv")

# load tral_perf panel for STRs that have STRs in their flanking regions where the consensus unit is the same or 
# part of the consensus unit (e.g. AT matches AAT)
df_neighbour_loose_match_panel <- load_str_panel("/Users/maxverbiest/PhD/data/str_panels/tral_and_perf_panel_flanks_cons_loose.tsv")

# add neighbour info to STR panel, construct neighbour_type column and drop columns for individual neighbour types
df_str_panel <- add_panel_info(df_str_panel, df_neighbour_panel, "has_neighbour")
df_str_panel <- add_panel_info(df_str_panel, df_neighbour_matching_panel, "has_matching_neighbour")
df_str_panel <- add_panel_info(df_str_panel, df_neighbour_loose_match_panel, "has_loose_matching_neighbour")
df_str_panel <- df_str_panel %>% 
  mutate(
    neighbour_type = case_when(
      has_matching_neighbour ~ "neighbour_match",
      has_loose_matching_neighbour ~ "neighbour_loose_match",
      has_neighbour ~ "neighbour",
      TRUE ~ "no_neighbour"
    )
  ) %>% 
  select(-c(has_matching_neighbour, has_loose_matching_neighbour, has_neighbour))

# write STR panel with meta information to file
write.table(df_str_panel, "/Users/maxverbiest/PhD/data/str_panels/tral_and_perf_panel_meta_info.tsv", sep="\t", quote=F, row.names=F)


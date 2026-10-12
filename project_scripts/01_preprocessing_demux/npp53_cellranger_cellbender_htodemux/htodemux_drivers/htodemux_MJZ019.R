source("/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/_htodemux_template.R")
hash_map <- list(
  "B301" = list(sample_id = "JZ201", model = "Npp53", treatment = "Castration + Enzalutamide", species = "mouse"),
  "B302" = list(sample_id = "JZ202", model = "Npp53", treatment = "Castration + Enzalutamide", species = "mouse"),
  "B303" = list(sample_id = "JZ203", model = "Npp53", treatment = "Castration + Enzalutamide", species = "mouse"),
  "B304" = list(sample_id = "JZ204", model = "Npp53", treatment = "Castration", species = "mouse"),
  "B306" = list(sample_id = "JZ205", model = "Npp53", treatment = "Castration", species = "mouse")
)
run_cohort("MJZ019", mt_pattern = "^mt-", hash_map = hash_map)

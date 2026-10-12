source("/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/_htodemux_template.R")
hash_map <- list(
  "B301" = list(sample_id = "JZ229", model = "Npp53", treatment = "Castration", species = "mouse"),
  "B302" = list(sample_id = "JZ230", model = "Npp53", treatment = "Castration", species = "mouse"),
  "B303" = list(sample_id = "JZ231", model = "Npp53", treatment = "Castration", species = "mouse"),
  "B304" = list(sample_id = "JZ232", model = "Npp53", treatment = "Castration", species = "mouse"),
  "B305" = list(sample_id = "JZ233", model = "Npp53", treatment = "Castration", species = "mouse"),
  "B306" = list(sample_id = "JZ234", model = "Npp53", treatment = "Castration", species = "mouse")
)
run_cohort("MJZ024", mt_pattern = "^mt-", hash_map = hash_map)

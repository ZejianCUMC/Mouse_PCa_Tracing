source("/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/_htodemux_template.R")
hash_map <- list(
  "B301" = list(sample_id = "JZ224", model = "Npp53", treatment = "Intact", species = "mouse"),
  "B302" = list(sample_id = "JZ225", model = "Npp53", treatment = "Intact", species = "mouse"),
  "B303" = list(sample_id = "JZ226", model = "Npp53", treatment = "Intact", species = "mouse"),
  "B304" = list(sample_id = "JZ227", model = "Npp53", treatment = "Intact", species = "mouse"),
  "B306" = list(sample_id = "JZ228", model = "Npp53", treatment = "Castration", species = "mouse")
)
run_cohort("MJZ023", mt_pattern = "^mt-", hash_map = hash_map)

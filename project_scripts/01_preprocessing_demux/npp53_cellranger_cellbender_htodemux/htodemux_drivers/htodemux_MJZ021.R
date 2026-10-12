source("/groups/ms3625_gp/zw2994/1-project/12-Jingbo/2-Tracing/scripts/_htodemux_template.R")
hash_map <- list(
  "B301" = list(sample_id = "JZ206", model = "Npp53", treatment = "Intact", species = "mouse"),
  "B302" = list(sample_id = "JZ207", model = "Npp53", treatment = "Intact", species = "mouse"),
  "B303" = list(sample_id = "JZ208", model = "Npp53", treatment = "Intact", species = "mouse"),
  "B304" = list(sample_id = "JZ209", model = "Npp53", treatment = "Intact", species = "mouse"),
  "B306" = list(sample_id = "JZ212", model = "Npp53", treatment = "Intact", species = "mouse")
)
run_cohort("MJZ021", mt_pattern = "^mt-", hash_map = hash_map)

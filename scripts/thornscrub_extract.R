library(terra)
library(sf)
library(tidyterra)
library(ggplot2)
library(rapr)

ecomap <- rast("E:/OneDrive - Texas A&M University - Kingsville/data/TX_Eco_Mapping/EcologicalMappingSystems_10mRaster.img")
ecomap
plot(ecomap)

# MCP for study areas?
mcps <- st_read("E:/GitProjects/GitprojectDroughtSelection/data/buffered_mcps.shp") %>% filter(Group == "area1")
mapview::mapview(mcps)
mcps <- st_transform(mcps, crs = crs(ecomap))

ecomap_stx <- crop(ecomap, mcps)
plot(ecomap_stx)
ecomap_stx
levels(ecomap_stx)
unique(values(ecomap_stx))

ecomap_stx[ecomap_stx == 0] <- NA

# Get the values actually present
used <- unique(values(ecomap_stx)[, 1])
used <- used[!is.na(used)]
# Keep only those levels
lev <- levels(ecomap_stx)[[1]]
lev <- lev[lev$value %in% used, ]
# Reassign the reduced levels table
levels(ecomap_stx) <- lev
freq(ecomap_stx)

# Plot stx eco map with mcps
ggplot() + geom_spatraster(data = ecomap_stx) +
  theme_minimal(base_size = 20) +
  geom_sf(data = mcps, fill = NA, linewidth = 1) +
  guides(fill = guide_legend(ncol = 2))

# Get categories
eco_cats <- cats(ecomap_stx)[[1]]

# Search for all woody canopy/brush types ocelots can utilize
ocelot_woody_pattern <- "Thornscrub|Thornforest|Shrubland|Ramadero|Woodland|Forest|Brush|Mesquite|Ebony|Cenizo"

# Drop non-woody land cover
exclude_non_woody <- "Grassland|Marsh|Prairie|Playa|Barren|Flat|Open Water"

target_legend <- eco_cats[
  (grepl(ocelot_woody_pattern, eco_cats$CommonName, ignore.case = TRUE) | 
   (eco_cats$value >= 7000 & eco_cats$value <= 7699) | 
   eco_cats$value %in% c(6806, 8306)) &
  !grepl(exclude_non_woody, eco_cats$CommonName, ignore.case = TRUE), 
]

target_values <- target_legend$value
ecomap_thorn <- ecomap_stx
ecomap_thorn[!(ecomap_thorn[] %in% target_values)] <- NA
ggplot() + geom_spatraster(data = ecomap_thorn) +
  theme_minimal(base_size = 20) +
  geom_sf(data = mcps, fill = NA, linewidth = 1) +
  guides(fill = guide_legend(ncol = 1))

ocelot_target_values <- ocelot_legend$value

# Build binary raster (1 = Woody/Thornscrub Cover, 0 = Other Land Cover)
thornscrub_binary <- ifel(is.na(ecomap_stx), NA, ifel(ecomap_stx %in% target_values, 1, 0)) %>% as.factor()

# 3. Plot binary map
ggplot() +
  geom_spatraster(data = thornscrub_binary) +
  scale_fill_manual(
    values = c("1" = "#2e6f40", "0" = "#d9c5b2"),
    labels = c("0" = "Other / Non-Woody", "1" = "Thornscrub & Woody Cover"),
    na.value = "transparent") +
  geom_sf(data = mcps, fill = NA, color = "black", linewidth = 0.8) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "right")

# Export as GeoTIFF with LZW compression to keep file size small
writeRaster(thornscrub_binary, "output/south_texas_thornscrub_binary.tif", 
  gdal = c("COMPRESS=LZW"), overwrite = TRUE)

# ============================================================================
# Compare to NLCD wood cover
# ============================================================================

nlcd <- rast("E:/OneDrive - Texas A&M University - Kingsville/data/nlcd/Annual_NLCD_LndCov_2022_CU_C1V1.tif")
plot(nlcd)

mcps <- st_transform(mcps, crs = crs(nlcd))
nlcd_stx <- crop(nlcd, mcps)
plot(nlcd_stx)
rm(nlcd)

# NLCD codes that represent woody structure
nlcd_woody_codes <- c(41, 42, 43, 52, 90)

# Binary NLCD raster (1 = Woody, 0 = Non-Woody)
nlcd_binary <- ifel(
  is.na(nlcd_stx),
  NA,
  ifel(nlcd_stx %in% nlcd_woody_codes, 1, 0)
) %>% as.factor()

# Plot NLCD Binary Map
ggplot() +
  geom_spatraster(data = nlcd_binary) +
  scale_fill_manual(
    values = c("1" = "#2e6f40", "0" = "#d9c5b2"),
    labels = c("0" = "Non-Woody", "1" = "Woody Cover (NLCD)"),
    na.value = "transparent",
    name = "NLCD Cover"
  ) +
  geom_sf(data = mcps, fill = NA, color = "black", linewidth = 0.8) +
  theme_minimal(base_size = 18)


nlcd_aligned <- project(nlcd_binary, thornscrub_binary, method = "near")
plot(nlcd_aligned)

diff_raster <- thornscrub_binary - nlcd_aligned

diff_raster <- as.numeric(thornscrub_binary) - as.numeric(nlcd_aligned)
diff_factor <- as.factor(diff_raster)

# 2. Plot simple subtraction map
ggplot() +
  geom_spatraster(data = diff_factor) +
  scale_fill_manual(
    values = c(
      "-1" = "#d9a752",  # NLCD Only
      "0"  = "#eaeaea",  # Agreement (Both present or both absent)
      "1"  = "#c0392b"   # Thornscrub Only
    ),
    labels = c(
      "-1" = "NLCD Only",
      "0"  = "Agreement (Combined)",
      "1"  = "Thornscrub Only"
    ),
    na.value = "transparent",
    name = "Difference"
  ) +
  geom_sf(data = mcps, fill = NA, color = "black", linewidth = 0.8) +
  theme_minimal(base_size = 18)

# Summarize agreement/disagreement pixel counts
diff_summary <- freq(diff_raster)
diff_summary$category <- c("NLCD only", "Agreement", "TxEco only")
diff_summary$area_ha <- (diff_summary$count * res(diff_raster)[1] * res(diff_raster)[2]) / 10000
diff_summary$area_km2 <- (diff_summary$count * res(diff_raster)[1] * res(diff_raster)[2]) * 1e-6

print(diff_summary[, c("category", "count", "area_ha", "area_km2")])

# ============================================================================
# Compare to RAP tree + shrub
# ============================================================================

# get_rap() needs crs 4326
mcp_4326 <- mcps %>% st_transform(crs = 4326) 
rap <- get_rap(mcp_4326, source = "rap-10m", product = "pft", years = 2024, verbose = TRUE)
rap

# Sum RAP tree and shrub cover bands
rap_total_woody <- rap$shrub + rap$tree

ggplot() +
  geom_spatraster(data = rap_total_woody) +
  geom_sf(data = mcps, fill = NA, color = "black", linewidth = 0.8) +
  theme_minimal(base_size = 18)


# Binary threshold (e.g., >= 30% total woody cover for inclusive ocelot cover)
rap_binary <- ifel(
  is.na(rap_total_woody),
  NA,
  ifel(rap_total_woody >= 20, 1, 0)
) %>% as.factor()

ggplot() +
  geom_spatraster(data = rap_binary) +
  scale_fill_manual(
    values = c("1" = "#2e6f40", "0" = "#d9c5b2"),
    labels = c("0" = "Non-Woody", "1" = "Woody Cover"),
    na.value = "transparent",
    name = "RAP Cover"
  ) +
  geom_sf(data = mcps, fill = NA, color = "black", linewidth = 0.8) +
  theme_minimal(base_size = 18)



























# # Get Names and IDs
# target_legend <- unique(target_legend)
# target_types  <- target_legend$VegID
# target_names  <- target_legend$CommonName

# message("--- POTENTIAL THORNSCRUB TYPES IN LEGEND ---")
# message("Total identified from classification system:", " ", length(target_types))

# # Calculate actual pixel counts across the entire tile
# all_freq <- freq(ecomap)

# # Merge categories with actual pixel counts
# target_summary <- merge(
#   target_legend[, c("VegID", "CommonName")], 
#   all_freq, 
#   by.x = "CommonName", 
#   by.y = "value", 
#   all.x = TRUE
# )

# target_summary$count[is.na(target_summary$count)] <- 0
# target_summary <- target_summary[order(target_summary$VegID), ]

# cat("--- AUDIT SUMMARY: ALL CANDIDATE TYPES ---\n")
# print(target_summary[, c("VegID", "CommonName", "count")], row.names = FALSE)

# # Get types with >0 pixels
# active_types <- target_summary[target_summary$count > 0, ]

# message("\n--- ACTIVE COVER TYPES IN THIS TILE ---\n")
# print(active_types[, c("VegID", "CommonName", "count")], row.names = FALSE)

# active_names <- active_types$CommonName

# # Get only target veg IDs and mask out everything else
# thornscrub_eco <- mask(ecomap, ecomap %in% target_summary$VegID, maskvalue = FALSE)

# # Plot to verify
# plot(target_eco)
# freq(target_eco)

# thornscrub_binary <- ifel(ecomap %in% target_summary$VegID, 1, NA)
# plot(thornscrub_binary)

# # Export as GeoTIFF with LZW compression to keep file size small
# writeRaster(thornscrub_binary, "output/south_texas_thornscrub_binary.tif", 
#   gdal = c("COMPRESS=LZW"), overwrite = TRUE)

# thornscrub_binary <- rast("output/south_texas_thornscrub_binary.tif")
# plot(thornscrub_binary)

# # thornscrub_specific_name <- c(7205, 7202, 7204, 7207, 7005, 7002, 7004, 6806, 8306)
# # target_legend <- eco_cats[eco_cats$value %in% thornscrub_specific_name,]
# # target_values <- target_legend$value

# # ecomap_thorn <- ecomap_stx
# # ecomap_thorn[!(ecomap_thorn[] %in% target_values)] <- NA
# # ggplot() + geom_spatraster(data = ecomap_thorn) +
# #   theme_minimal(base_size = 20) +
# #   geom_sf(data = mcps, fill = NA, linewidth = 1) +
# #   guides(fill = guide_legend(ncol = 1))

# # # Also look for associated names and vegetation
# # thorn_pattern <- "Thornscrub|Thornforest|Shrubland|Ramadero|Blackbrush|Ebony|Cenizo|Granjeno|Tamaulipan"
# # target_legend <- eco_cats[grepl(thorn_pattern, eco_cats$CommonName, ignore.case = TRUE),]
# target_values <- target_legend$value

# ecomap_thorn <- ecomap_stx
# ecomap_thorn[!(ecomap_thorn[] %in% target_values)] <- NA
# ggplot() + geom_spatraster(data = ecomap_thorn) +
#   theme_minimal(base_size = 20) +
#   geom_sf(data = mcps, fill = NA, linewidth = 1) +
#   guides(fill = guide_legend(ncol = 1))

# target_legend <- eco_cats[grepl(thorn_pattern, eco_cats$CommonName, ignore.case = TRUE) | 
#   (eco_cats$value >= 7000 & eco_cats$value <= 7699) | # Tamaulipan / S. TX thornscrub IDs
#   eco_cats$value %in% c(6806, 8306),                 # Salty Thornscrub & Trans-Pecos Shrubland
#   ]
# target_values <- target_legend$value

# ecomap_thorn2 <- ecomap_stx
# ecomap_thorn2[!(ecomap_thorn2[] %in% target_values)] <- NA
# ggplot() + geom_spatraster(data = ecomap_thorn2) +
#   theme_minimal(base_size = 20) +
#   geom_sf(data = mcps, fill = NA, linewidth = 1) +
#   guides(fill = guide_legend(ncol = 1))

# # Binary raster
# ecomap_binary <- ifel(
#   is.na(ecomap_stx),
#   NA,
#   ifel(ecomap_stx %in% target_values, 1, 0)
# ) %>% as.factor()

# ggplot() + geom_spatraster(data = ecomap_binary) +
#   scale_fill_manual(
#     values = c("1" = "lightblue", "0" = "darkblue"),
#     na.value = "transparent"
#   ) +  theme_minimal(base_size = 20) +
#   #geom_sf(data = mcps, fill = NA, linewidth = 1) +
#   guides(fill = guide_legend(ncol = 1))
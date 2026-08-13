library(terra)

ecomap <- rast("E:/OneDrive - Texas A&M University - Kingsville/data/TX_Eco_Mapping/South_Texas_Gulf_Eco_Mapping.tif")
ecomap
plot(ecomap)

# Get categories
eco_cats <- cats(ecomap)[[1]]

# Also look for associated names and vegetation
thorn_pattern <- "Thornscrub|Thornforest|Shrubland|Ramadero|Blackbrush|Ebony|Cenizo|Granjeno|Tamaulipan"

target_legend <- eco_cats[
  grepl(thorn_pattern, eco_cats$CommonName, ignore.case = TRUE) | 
  (eco_cats$VegID >= 7000 & eco_cats$VegID <= 7699) | # Tamaulipan / S. TX thornscrub IDs
  eco_cats$VegID %in% c(6806, 8306),                 # Salty Thornscrub & Trans-Pecos Shrubland
]

# Get Names and IDs
target_legend <- unique(target_legend)
target_types  <- target_legend$VegID
target_names  <- target_legend$CommonName

message("--- POTENTIAL THORNSCRUB TYPES IN LEGEND ---")
message("Total identified from classification system:", " ", length(target_types))

# Calculate actual pixel counts across the entire tile
all_freq <- freq(ecomap)

# Merge categories with actual pixel counts
target_summary <- merge(
  target_legend[, c("VegID", "CommonName")], 
  all_freq, 
  by.x = "CommonName", 
  by.y = "value", 
  all.x = TRUE
)

target_summary$count[is.na(target_summary$count)] <- 0
target_summary <- target_summary[order(target_summary$VegID), ]

cat("--- AUDIT SUMMARY: ALL CANDIDATE TYPES ---\n")
print(target_summary[, c("VegID", "CommonName", "count")], row.names = FALSE)

# Get types with >0 pixels
active_types <- target_summary[target_summary$count > 0, ]

message("\n--- ACTIVE COVER TYPES IN THIS TILE ---\n")
print(active_types[, c("VegID", "CommonName", "count")], row.names = FALSE)

active_names <- active_types$CommonName

# Get only target veg IDs and mask out everything else
thornscrub_eco <- mask(ecomap, ecomap %in% target_summary$VegID, maskvalue = FALSE)

# Plot to verify
plot(target_eco)
freq(target_eco)

thornscrub_binary <- ifel(ecomap %in% target_summary$VegID, 1, NA)
plot(thornscrub_binary)

# Export as GeoTIFF with LZW compression to keep file size small
writeRaster(thornscrub_binary, "output/south_texas_thornscrub_binary.tif", 
  gdal = c("COMPRESS=LZW"), overwrite = TRUE)


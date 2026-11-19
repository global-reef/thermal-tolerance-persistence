install.packages(c("sf", "ggplot2", "dplyr", "osmdata", "ggrepel", "ggspatial", "patchwork"))

### 01. SETUP ####

library(sf)
library(ggplot2)
library(dplyr)
library(osmdata)
library(ggrepel)
library(ggspatial)
library(patchwork)

theme_set(theme_void())

### 02. SITE COORDINATES ####

sites <- tibble::tribble(
  ~site,              ~lat,        ~lon,
  "Shallow donor",    10.0637439,  99.8229428,
  "Deep donor",       10.0639909,  99.8243499,
  "Transplant site",  10.0629401,  99.8241716
)

sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)

### 03. MAP EXTENTS ####

# Close-up around Chalok Bay
bbox_main <- c(
  xmin = 99.8205,
  ymin = 10.0618,
  xmax = 99.8268,
  ymax = 10.0649
)

# Island-scale locator
bbox_island <- c(
  xmin = 99.790,
  ymin = 10.020,
  xmax = 99.870,
  ymax = 10.115
)

main_box_sf <- st_as_sfc(st_bbox(bbox_main, crs = 4326))
island_box_sf <- st_as_sfc(st_bbox(bbox_island, crs = 4326))

### 04. RETRIEVE KOH TAO POLYGON FROM OSM ####

# Query island polygons
island_query <- opq(bbox = bbox_island, timeout = 60) %>%
  add_osm_feature(key = "place", value = "island") %>%
  osmdata_sf()

island_polys <- bind_rows(
  island_query$osm_polygons,
  island_query$osm_multipolygons
) %>%
  st_make_valid()

# Keep Koh Tao / Ko Tao only
koh_tao <- island_polys %>%
  filter(name %in% c("Ko Tao", "Koh Tao"))

# If multiple features returned, keep the largest polygon
if (nrow(koh_tao) > 1) {
  koh_tao <- koh_tao %>%
    mutate(area = st_area(geometry)) %>%
    slice_max(area, n = 1) %>%
    select(-area)
}

# Crop for close-up panel
koh_tao_main <- st_crop(koh_tao, st_bbox(main_box_sf))

### 05. LABEL POSITIONS ####

site_labels <- sites_sf %>%
  mutate(
    label_x = c(99.8218, 99.8252, 99.8251),
    label_y = c(10.0644, 10.0645, 10.0623)
  ) %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry()

### 06. MAIN PANEL: CHALOK BAY ####

p_main <- ggplot() +
  geom_sf(data = koh_tao_main, fill = "grey90", colour = "black", linewidth = 0.35) +
  geom_sf(
    data = sites_sf,
    aes(shape = site, fill = site),
    colour = "black",
    size = 3,
    stroke = 0.5,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = site_labels,
    aes(X, Y, label = site),
    size = 3.2,
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.3,
    point.padding = 0.25,
    segment.size = 0.3
  ) +
  scale_shape_manual(values = c(
    "Shallow donor" = 21,
    "Deep donor" = 24,
    "Transplant site" = 22
  )) +
  scale_fill_manual(values = c(
    "Shallow donor" = "white",
    "Deep donor" = "grey40",
    "Transplant site" = "black"
  )) +
  coord_sf(
    xlim = c(bbox_main["xmin"], bbox_main["xmax"]),
    ylim = c(bbox_main["ymin"], bbox_main["ymax"]),
    expand = FALSE
  ) +
  annotation_scale(location = "br", width_hint = 0.28, text_cex = 0.7) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_orienteering(
      fill = c("white", "black"),
      line_col = "black"
    ),
    height = unit(0.8, "cm"),
    width = unit(0.8, "cm")
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(5, 5, 5, 5)
  )

### 07. LOCATOR PANEL: KOH TAO ####

study_point <- sites_sf %>%
  filter(site == "Transplant site")

p_inset <- ggplot() +
  geom_sf(data = koh_tao, fill = "grey90", colour = "black", linewidth = 0.35) +
  geom_sf(data = main_box_sf, fill = NA, colour = "black", linewidth = 0.45, linetype = "dashed") +
  geom_sf(data = study_point, shape = 21, fill = "black", colour = "black", size = 2) +
  annotate("text", x = 99.833, y = 10.087, label = "Koh Tao", size = 3.4) +
  coord_sf(
    xlim = c(bbox_island["xmin"], bbox_island["xmax"]),
    ylim = c(bbox_island["ymin"], bbox_island["ymax"]),
    expand = FALSE
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(2, 2, 2, 2)
  )

### 08. COMBINE ####

final_map <- p_main +
  inset_element(
    p_inset,
    left = 0.03, bottom = 0.60, right = 0.33, top = 0.98,
    align_to = "panel"
  )

print(final_map)

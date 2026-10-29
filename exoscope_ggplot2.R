# Author: Kirti Anand Nalage
# Project: ExoScope: Exoplanet Discovery Trends and Earth-Similarity Analysis using R
# ------------------------------------------------------------------
# Data: NASA Exoplanet Archive, Planetary Systems Composite Parameters
# table (pscomppars) -- one row per confirmed planet.
#
# How to run: keep this script in an empty folder, open it in RStudio,
# click Source. First run downloads the data and saves a local snapshot
# (exoplanets_raw.csv). All charts go to the plots/ folder.
# ------------------------------------------------------------------

library(dplyr)
library(ggplot2)

dir.create("plots", showWarnings = FALSE)
theme_set(theme_minimal(base_size = 12))

# ---- 1. Load data (local snapshot first) ----
snapshot <- "exoplanets_raw.csv"

if (file.exists(snapshot)) {
  planets_raw <- read.csv(snapshot, stringsAsFactors = FALSE)
} else {
  url <- paste0(
    "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+",
    "pl_name,hostname,discoverymethod,disc_year,pl_orbper,pl_rade,",
    "pl_bmasse,pl_eqt,st_teff,sy_dist+from+pscomppars&format=csv"
  )
  planets_raw <- read.csv(url, stringsAsFactors = FALSE)
  write.csv(planets_raw, snapshot, row.names = FALSE)
  cat("Data downloaded on:", format(Sys.Date()), "(write this date in the report)\n")
}

needed <- c("pl_name", "hostname", "discoverymethod", "disc_year", "pl_orbper",
            "pl_rade", "pl_bmasse", "pl_eqt", "st_teff", "sy_dist")
missing_cols <- setdiff(needed, names(planets_raw))
if (length(missing_cols) > 0) {
  stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
}

# ---- 2. Clean ----
# Methods with fewer than 30 planets are grouped as "Other" so charts stay readable
planets <- planets_raw %>%
  select(all_of(needed)) %>%
  filter(!is.na(disc_year)) %>%
  add_count(discoverymethod, name = "n_method") %>%
  mutate(method_group = ifelse(n_method < 30, "Other", discoverymethod))

cat("\nTotal planets:", nrow(planets), "\n")
cat("Missing values per column:\n")
print(colSums(is.na(planets[, needed])))

# ---- 3. Discovery trend ----
trend <- planets %>% count(disc_year, method_group)

p1 <- ggplot(trend, aes(x = disc_year, y = n, fill = method_group)) +
  geom_col() +
  scale_fill_viridis_d(option = "D", name = "Discovery method") +
  labs(title = "Exoplanet discoveries per year, by method",
       x = "Discovery year", y = "Planets discovered")
ggsave("plots/01_discoveries_per_year.png", p1, width = 9, height = 5.5, dpi = 150)

cum <- planets %>%
  count(disc_year) %>%
  arrange(disc_year) %>%
  mutate(total = cumsum(n))

p2 <- ggplot(cum, aes(x = disc_year, y = total)) +
  geom_line(color = "#264653") +
  geom_point(color = "#e76f51") +
  labs(title = "Confirmed exoplanets, cumulative count",
       x = "Discovery year", y = "Total confirmed planets")
ggsave("plots/02_cumulative_discoveries.png", p2, width = 9, height = 5.5, dpi = 150)

# ---- 4. Size vs orbit map ----
sz <- planets %>%
  filter(!is.na(pl_rade), !is.na(pl_orbper), pl_rade > 0, pl_orbper > 0)

cat("\nPlanets with both radius and orbital period:", nrow(sz), "\n")

p3 <- ggplot(sz, aes(x = pl_orbper, y = pl_rade, color = method_group)) +
  geom_point(alpha = 0.6, size = 1.4) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_d(option = "D", name = "Discovery method") +
  annotate("point", x = 365.25, y = 1, shape = 8, size = 4) +
  annotate("text", x = 365.25, y = 0.8, label = "Earth", size = 3.5) +
  annotate("point", x = 4333, y = 11.2, shape = 8, size = 4) +
  annotate("text", x = 4333, y = 9, label = "Jupiter", size = 3.5) +
  labs(title = "Planet size vs orbital period (log scales)",
       x = "Orbital period (days)", y = "Radius (Earth radii)")
ggsave("plots/03_size_vs_orbit.png", p3, width = 9, height = 6, dpi = 150)

# ---- 5. Clustering (k-means on log radius + log period) ----
cl_data <- sz %>%
  mutate(log_rade = log10(pl_rade), log_per = log10(pl_orbper))

X <- scale(cl_data[, c("log_rade", "log_per")])

set.seed(42)
elbow <- data.frame(
  k = 1:8,
  wss = sapply(1:8, function(k) kmeans(X, centers = k, nstart = 25)$tot.withinss)
)

p4 <- ggplot(elbow, aes(x = k, y = wss)) +
  geom_line(color = "#264653") +
  geom_point(size = 3, color = "#e76f51") +
  scale_x_continuous(breaks = 1:8) +
  labs(title = "Elbow plot: choosing the number of clusters",
       x = "Number of clusters (k)", y = "Total within-cluster sum of squares")
ggsave("plots/04_elbow_plot.png", p4, width = 7, height = 5, dpi = 150)

K <- 4   # change after looking at the elbow plot

km <- kmeans(X, centers = K, nstart = 25)
cl_data$cluster_raw <- km$cluster

# Order clusters by median radius so labels are stable: C1 = smallest planets
ord <- cl_data %>%
  group_by(cluster_raw) %>%
  summarise(med = median(pl_rade), .groups = "drop") %>%
  arrange(med) %>%
  mutate(cluster_rank = row_number())

cl_data <- cl_data %>% left_join(select(ord, cluster_raw, cluster_rank), by = "cluster_raw")

# Rename these after reading cluster_summary below (e.g. "Rocky", "Gas giants")
cluster_names <- c("C1 (smallest)", paste0("C", 2:(K - 1)), paste0("C", K, " (largest)"))
cl_data$cluster <- factor(cl_data$cluster_rank, levels = 1:K, labels = cluster_names)

cluster_summary <- cl_data %>%
  group_by(cluster) %>%
  summarise(n = n(),
            median_radius_earth = round(median(pl_rade), 2),
            median_period_days = round(median(pl_orbper), 2),
            .groups = "drop")

cat("\nCluster summary:\n")
print(cluster_summary)
write.csv(cluster_summary, "cluster_summary.csv", row.names = FALSE)

p5 <- ggplot(cl_data, aes(x = pl_orbper, y = pl_rade, color = cluster)) +
  geom_point(alpha = 0.6, size = 1.4) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_d(option = "C", name = "Cluster") +
  labs(title = paste0("K-means clusters of exoplanets (k = ", K, ")"),
       x = "Orbital period (days)", y = "Radius (Earth radii)")
ggsave("plots/05_clusters.png", p5, width = 9, height = 6, dpi = 150)

# ---- 6. Earth-similarity score (ESI-lite) ----
# Simplified version of the Earth Similarity Index using only the two
# parameters available in our table: radius and equilibrium temperature.
# Earth: radius = 1 Earth radius, equilibrium temperature ~ 255 K.
# This is a screening score only. It does NOT say a planet is habitable.
earth_r <- 1
earth_t <- 255
sim <- function(x, x0) 1 - abs(x - x0) / (x + x0)

esi <- planets %>%
  filter(!is.na(pl_rade), !is.na(pl_eqt), pl_rade > 0, pl_eqt > 0) %>%
  mutate(sim_radius = sim(pl_rade, earth_r),
         sim_temp = sim(pl_eqt, earth_t),
         esi_lite = sim_radius^(0.57 / 2) * sim_temp^(5.58 / 2))

cat("\nPlanets with both radius and equilibrium temperature:", nrow(esi), "\n")

top10 <- esi %>%
  arrange(desc(esi_lite)) %>%
  slice_head(n = 10) %>%
  select(pl_name, hostname, pl_rade, pl_eqt, esi_lite) %>%
  mutate(esi_lite = round(esi_lite, 3))

cat("\nTop 10 most Earth-like (ESI-lite):\n")
print(as.data.frame(top10))
write.csv(top10, "top10_earth_similar.csv", row.names = FALSE)

p6 <- ggplot(top10, aes(x = reorder(pl_name, esi_lite), y = esi_lite)) +
  geom_col(fill = "#2a9d8f") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Top 10 planets by Earth-similarity (ESI-lite)",
       x = NULL, y = "ESI-lite score (1 = identical to Earth on these two traits)")
ggsave("plots/06_top10_earth_similar.png", p6, width = 8, height = 5.5, dpi = 150)

cat("\nDone. Charts are in the plots/ folder.\n")

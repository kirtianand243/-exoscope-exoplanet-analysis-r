# Author: Kirti Anand Nalage
# Project: ExoScope: Exoplanet Discovery Trends and Earth-Similarity Analysis using R
# ------------------------------------------------------------------
# BASE R VERSION: needs NO extra packages (no dplyr, no ggplot2).
# Data: NASA Exoplanet Archive, Planetary Systems Composite Parameters
# (pscomppars) -- one row per confirmed planet.
#
# How to run:
#  1. Put this script and exoplanets_raw.csv in the SAME folder
#  2. RStudio: Session > Set Working Directory > To Source File Location
#  3. Click Source
# Charts are saved in the plots/ folder.
# ------------------------------------------------------------------

cat("Working folder:", getwd(), "\n")
dir.create("plots", showWarnings = FALSE)

# helper: save a chart to file AND show it on screen
make_plot <- function(file, draw, w = 1200, h = 700) {
  png(file, width = w, height = h, res = 130)
  draw()
  invisible(dev.off())
  draw()
}

# ---- 1. Load data (local file first) ----
snapshot <- "exoplanets_raw.csv"

if (file.exists(snapshot)) {
  raw <- read.csv(snapshot, stringsAsFactors = FALSE)
} else {
  url <- paste0(
    "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+",
    "pl_name,hostname,discoverymethod,disc_year,pl_orbper,pl_rade,",
    "pl_bmasse,pl_eqt,st_teff,sy_dist+from+pscomppars&format=csv"
  )
  raw <- tryCatch(read.csv(url, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(raw)) {
    stop("Download failed. Open the link in a browser, save the CSV as ",
         "exoplanets_raw.csv inside this folder: ", getwd())
  }
  write.csv(raw, snapshot, row.names = FALSE)
  cat("Data downloaded on:", format(Sys.Date()), "(write this date in the report)\n")
}

needed <- c("pl_name", "hostname", "discoverymethod", "disc_year", "pl_orbper",
            "pl_rade", "pl_bmasse", "pl_eqt", "st_teff", "sy_dist")
missing_cols <- setdiff(needed, names(raw))
if (length(missing_cols) > 0) {
  stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
}

# ---- 2. Clean ----
p <- raw[, needed]
p <- p[!is.na(p$disc_year), ]

# Methods with fewer than 30 planets are grouped as "Other"
method_counts <- table(p$discoverymethod)
small_methods <- names(method_counts)[method_counts < 30]
p$method_group <- ifelse(p$discoverymethod %in% small_methods, "Other", p$discoverymethod)

cat("\nTotal planets:", nrow(p), "\n")
cat("Missing values per column:\n")
print(colSums(is.na(p[, needed])))

# ---- 3. Discovery trend ----
tt <- table(p$method_group, p$disc_year)
method_cols <- rainbow(nrow(tt), s = 0.7, v = 0.85)

make_plot("plots/01_discoveries_per_year.png", function() {
  barplot(tt, col = method_cols, las = 2, cex.names = 0.8,
          main = "Exoplanet discoveries per year, by method",
          xlab = "Discovery year", ylab = "Planets discovered",
          legend.text = rownames(tt),
          args.legend = list(x = "topleft", cex = 0.75, bty = "n"))
})

per_year <- table(p$disc_year)
years <- as.numeric(names(per_year))

make_plot("plots/02_cumulative_discoveries.png", function() {
  plot(years, cumsum(per_year), type = "o", pch = 16, col = "#264653",
       main = "Confirmed exoplanets, cumulative count",
       xlab = "Discovery year", ylab = "Total confirmed planets")
})

# ---- 4. Size vs orbit map ----
ok <- !is.na(p$pl_rade) & !is.na(p$pl_orbper) & p$pl_rade > 0 & p$pl_orbper > 0
sz <- p[ok, ]
cat("\nPlanets with both radius and orbital period:", nrow(sz), "\n")

groups <- sort(unique(sz$method_group))
group_cols <- rainbow(length(groups), s = 0.7, v = 0.85)
pt_col <- adjustcolor(group_cols[match(sz$method_group, groups)], alpha.f = 0.5)

make_plot("plots/03_size_vs_orbit.png", function() {
  plot(sz$pl_orbper, sz$pl_rade, log = "xy", pch = 16, col = pt_col, cex = 0.8,
       xlim = range(c(sz$pl_orbper, 365.25, 4333)),
       main = "Planet size vs orbital period (log scales)",
       xlab = "Orbital period (days)", ylab = "Radius (Earth radii)")
  points(365.25, 1, pch = 8, cex = 2)
  text(365.25, 0.75, "Earth")
  points(4333, 11.2, pch = 8, cex = 2)
  text(4333, 8.5, "Jupiter")
  legend("topright", legend = groups, col = group_cols, pch = 16, cex = 0.7, bty = "n")
})

# ---- 5. Clustering (k-means on log radius + log period) ----
X <- scale(cbind(log10(sz$pl_rade), log10(sz$pl_orbper)))

set.seed(42)
wss <- sapply(1:8, function(k) kmeans(X, centers = k, nstart = 25)$tot.withinss)

make_plot("plots/04_elbow_plot.png", function() {
  plot(1:8, wss, type = "o", pch = 16, col = "#264653",
       main = "Elbow plot: choosing the number of clusters",
       xlab = "Number of clusters (k)", ylab = "Total within-cluster sum of squares")
}, w = 900, h = 650)

K <- 4   # change after looking at the elbow plot (use 2 or more)

km <- kmeans(X, centers = K, nstart = 25)

# Order clusters by median radius so labels are stable: C1 = smallest planets
med_radius <- tapply(sz$pl_rade, km$cluster, median)
rank_of <- rank(med_radius, ties.method = "first")
sz$cluster_rank <- as.integer(rank_of[km$cluster])

# Rename these after reading cluster_summary below (e.g. "Rocky", "Gas giants")
cluster_names <- paste0("C", 1:K)
cluster_names[1] <- paste(cluster_names[1], "(smallest)")
cluster_names[K] <- paste(cluster_names[K], "(largest)")

cluster_summary <- data.frame(
  cluster = cluster_names,
  n = as.vector(table(sz$cluster_rank)),
  median_radius_earth = round(as.vector(tapply(sz$pl_rade, sz$cluster_rank, median)), 2),
  median_period_days = round(as.vector(tapply(sz$pl_orbper, sz$cluster_rank, median)), 2)
)

cat("\nCluster summary:\n")
print(cluster_summary)
write.csv(cluster_summary, "cluster_summary.csv", row.names = FALSE)

cl_cols <- rainbow(K, s = 0.8, v = 0.8)

make_plot("plots/05_clusters.png", function() {
  plot(sz$pl_orbper, sz$pl_rade, log = "xy", pch = 16, cex = 0.8,
       col = adjustcolor(cl_cols[sz$cluster_rank], alpha.f = 0.5),
       main = paste0("K-means clusters of exoplanets (k = ", K, ")"),
       xlab = "Orbital period (days)", ylab = "Radius (Earth radii)")
  legend("topright", legend = cluster_names, col = cl_cols, pch = 16, cex = 0.8, bty = "n")
})

# ---- 6. Earth-similarity score (ESI-lite) ----
# Simplified version of the Earth Similarity Index using only the two
# parameters available in our table: radius and equilibrium temperature.
# Earth: radius = 1 Earth radius, equilibrium temperature ~ 255 K.
# This is a screening score only. It does NOT say a planet is habitable.
sim <- function(x, x0) 1 - abs(x - x0) / (x + x0)

ok2 <- !is.na(p$pl_rade) & !is.na(p$pl_eqt) & p$pl_rade > 0 & p$pl_eqt > 0
e <- p[ok2, ]
e$esi_lite <- sim(e$pl_rade, 1)^(0.57 / 2) * sim(e$pl_eqt, 255)^(5.58 / 2)

cat("\nPlanets with both radius and equilibrium temperature:", nrow(e), "\n")

e <- e[order(-e$esi_lite), ]
top10 <- head(e[, c("pl_name", "hostname", "pl_rade", "pl_eqt", "esi_lite")], 10)
top10$esi_lite <- round(top10$esi_lite, 3)

cat("\nTop 10 most Earth-like (ESI-lite):\n")
print(top10, row.names = FALSE)
write.csv(top10, "top10_earth_similar.csv", row.names = FALSE)

make_plot("plots/06_top10_earth_similar.png", function() {
  op <- par(mar = c(5, 11, 4, 2))
  barplot(rev(top10$esi_lite), names.arg = rev(top10$pl_name), horiz = TRUE,
          las = 1, xlim = c(0, 1), col = "#2a9d8f",
          main = "Top 10 planets by Earth-similarity (ESI-lite)",
          xlab = "ESI-lite score (1 = identical to Earth on these two traits)")
  par(op)
}, w = 1100, h = 700)

cat("\nDone. Charts are in the plots/ folder.\n")

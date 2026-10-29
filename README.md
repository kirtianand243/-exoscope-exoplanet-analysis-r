# ExoScope: Exoplanet Discovery Trends and Earth-Similarity Analysis using R

An R data-analysis project on NASA's public table of confirmed exoplanets. It cleans the data, charts how discoveries have grown, groups planets into types with k-means clustering, and ranks known planets by how similar they are to Earth on two measurable traits.

Built as a micro project for the course *Data Science using R* at Annasaheb Dange College of Engineering and Technology, Ashta.

**Author:** Kirti Anand Nalage (B.Tech, Artificial Intelligence and Data Science)

## What the project answers

1. **How have exoplanet discoveries grown, and by which methods?** Cumulative and yearly counts, split by detection method.
2. **What kinds of planets are known?** Planet size vs orbital period, and k-means clusters (k = 4) on log radius and log period.
3. **Which known planets look most like Earth?** A simplified Earth Similarity Index (ESI-lite) using radius and equilibrium temperature.

## Results

Data downloaded in October 2026: about 6,450 confirmed planets.

| Discoveries over time | Size vs orbital period |
|---|---|
| ![Cumulative discoveries](figures/1_cumulative_discoveries.png) | ![Size vs orbit](figures/2_size_vs_orbit.png) |

| K-means clusters | Top 10 by ESI-lite |
|---|---|
| ![Clusters](figures/3_kmeans_clusters.png) | ![Top 10](figures/4_top10_earth_similar.png) |

- Discoveries grew slowly until about 2013, then jumped around 2014 and 2016 (years commonly linked to large batches of planets from NASA's Kepler mission).
- The known planets are strongly shaped by detection methods. Large planets on short orbits are easy to find, so the region around an Earth-like planet (Earth-sized, one-year orbit) is sparse in the chart.
- The four clusters can be read as small rocky / super-Earth planets, sub-Neptune and Neptune-sized planets, hot giants (hot-Jupiter type) and cool long-period giants.
- Top ESI-lite candidates, in order: Wolf 1069 b, TRAPPIST-1 e, TOI-700 d, GJ 1061 d, TOI-700 e, Kepler-442 b, Teegarden's Star b, Kepler-1649 c, Kepler-452 b, GJ 3378 b.

## How the Earth-similarity score works

For a property x with Earth value x0, similarity = 1 - |x - x0| / (x + x0). ESI-lite multiplies the two similarities, each raised to (published ESI weight / 2):

```
ESI-lite = (radius similarity)^(0.57 / 2) * (temperature similarity)^(5.58 / 2)
```

Earth is taken as radius = 1 Earth radius and equilibrium temperature = 255 K. A score of 1 means identical to Earth on these two traits. **This is a screening score only. It does not claim that any planet is habitable.**

## Files

| File | What it is |
|---|---|
| `exoscope.R` | Main script. Base R only, no extra packages needed. |
| `exoscope_ggplot2.R` | Alternative version of the same analysis using `dplyr` and `ggplot2`. |
| `figures/` | Result charts shown above. |

## How to run

1. Install R (and optionally RStudio).
2. Download this repository and keep the script in an empty folder.
3. In RStudio: open `exoscope.R`, then **Session > Set Working Directory > To Source File Location**, then click **Source**.
4. On the first run the script downloads the data from the NASA Exoplanet Archive and saves it as `exoplanets_raw.csv`. If the download fails, open the URL inside the script in a browser, save the CSV as `exoplanets_raw.csv` in the same folder, and run again.
5. Charts are saved in a `plots/` folder, and two result tables (`cluster_summary.csv`, `top10_earth_similar.csv`) are written next to the script.

For `exoscope_ggplot2.R`, install the packages once first: `install.packages(c("dplyr", "ggplot2"))`.

## Limitations

- Radius, mass and temperature are missing for many planets, so each analysis uses only the rows that have the variables it needs. These rows are not a random sample of all planets.
- The archive shows what telescopes can detect, not what exists (detection bias).
- K-means always returns k clusters, even when the data is continuous. The choice of k = 4 was guided by the elbow method and by interpretability.
- ESI-lite uses only two properties. Equilibrium temperature is a theoretical value, and radius alone does not tell whether a planet is rocky.
- The archive is updated regularly, so numbers depend on the download date.

## Data source and acknowledgement

This project has made use of the NASA Exoplanet Archive, which is operated by the California Institute of Technology, under contract with the National Aeronautics and Space Administration under the Exoplanet Exploration Program. Table used: Planetary Systems Composite Parameters (`pscomppars`). https://exoplanetarchive.ipac.caltech.edu

## References

- R Core Team, *R: A Language and Environment for Statistical Computing*. https://www.R-project.org/
- Schulze-Makuch, D. et al., "A Two-Tiered Approach to Assessing the Habitability of Exoplanets," *Astrobiology*, 11(10), 2011.
- Hartigan, J. A. and Wong, M. A., "Algorithm AS 136: A K-Means Clustering Algorithm," *Applied Statistics*, 28(1), 1979.

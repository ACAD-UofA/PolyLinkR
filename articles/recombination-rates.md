# Incorporating Recombination Rates into Polygenic Linkage Analysis

``` r
library(polylinkR)
library(data.table)
```

## 1. What Recombination Rates Add to the Analysis

Recombination rates are fundamental to understanding genetic linkage and
linkage disequilibrium (LD) patterns. Incorporating them into polygenic
pathway enrichment analysis provides several advantages:

### 1.1 Biological Rationale

**Recombination shapes LD patterns:** - Hotspots of recombination break
up LD blocks - Cold regions maintain long-range LD - The relationship
between physical distance (bp) and genetic distance (cM) varies across
the genome

**Why physical distance alone is insufficient:** - 1 Mb in a
recombination hotspot ≠ 1 Mb in a cold region (genetically) - Genes
500kb apart in a hotspot may be less correlated than genes 2Mb apart in
a cold region - Accounting for recombination improves estimates of
genetic autocorrelation

### 1.2 Impact on Analysis

| Without Recombination Rates                     | With Recombination Rates                                    |
|-------------------------------------------------|-------------------------------------------------------------|
| Physical distance (bp) used for autocorrelation | Genetic distance (cM) used for autocorrelation              |
| May underestimate LD in cold regions            | More accurate LD estimation across genome                   |
| May overestimate independence in hotspots       | Corrects for recombination breakpoints                      |
| Uniform kernel weights across chromosomes       | Distance-adjusted weights reflecting true genetic proximity |

### 1.3 Key Applications

1.  **Autocorrelation estimation**: Distance between genes measured in
    centiMorgans (cM) instead of base pairs
2.  **Gene set rescaling**: Decorrelation weights based on genetic
    rather than physical distance
3.  **Crossover-aware null models**: Permutation nulls that respect
    recombination boundaries

## 2. Loading Recombination Map Data (HapMap Format, Averaged from deCODE)

polylinkR includes recombination rate data from Bherer et al. (2017),
which provides sex-averaged recombination rates derived from deCODE
genetics data.

### 2.1 Included Data

``` r
# Load the included recombination rate data
data(human_sexavg_RR_Bherer2017)

# Check the structure
str(human_sexavg_RR_Bherer2017)
head(human_sexavg_RR_Bherer2017)
```

### 2.2 Data Format

The recombination rate file follows HapMap format with these columns:

| Column | Description                 | Units                      |
|--------|-----------------------------|----------------------------|
| `chr`  | Chromosome number           | Integer (1-22, X=23, Y=24) |
| `pos`  | Physical position           | Base pairs (bp)            |
| `rate` | Recombination rate          | cM per Mb                  |
| `cM`   | Cumulative genetic distance | CentiMorgans (cM)          |

### 2.3 Converting Physical to Genetic Distance

The `cM` column provides cumulative genetic map position. You can
calculate genetic distances between positions:

``` r
# Example: Calculate genetic distance between two positions on chromosome 1
rr_chr1 <- human_sexavg_RR_Bherer2017[chr == 1]

# Find genetic positions for two physical positions
pos1 <- 1000000
pos2 <- 5000000

# Interpolate to get genetic distances
cM1 <- approx(rr_chr1$pos, rr_chr1$cM, xout = pos1)$y
cM2 <- approx(rr_chr1$pos, rr_chr1$cM, xout = pos2)$y

if (!is.na(cM1) && !is.na(cM2)) {
  genetic_distance <- abs(cM2 - cM1)
  physical_distance <- abs(pos2 - pos1) / 1e6  # in Mb
  
  cat("Physical distance:", physical_distance, "Mb\n")
  cat("Genetic distance:", round(genetic_distance, 4), "cM\n")
  cat("Recombination rate in this region:", 
      round(genetic_distance / physical_distance, 4), "cM/Mb\n")
}
```

### 2.4 Visualizing Recombination Landscapes

``` r
# Plot recombination rate across chromosome 1
rr_chr1_plot <- human_sexavg_RR_Bherer2017[chr == 1]

# Sample every 100th point for faster plotting
plot_indices <- seq(1, nrow(rr_chr1_plot), by = 100)

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# Recombination rate
plot(rr_chr1_plot$pos[plot_indices] / 1e6, 
     rr_chr1_plot$rate[plot_indices],
     type = "l", col = "blue", lwd = 0.5,
     xlab = "Position (Mb)", ylab = "cM/Mb",
     main = "Recombination Rate across Chromosome 1")
abline(h = mean(rr_chr1_plot$rate, na.rm = TRUE), col = "red", lty = 2)

# Cumulative genetic map
plot(rr_chr1_plot$pos[plot_indices] / 1e6, 
     rr_chr1_plot$cM[plot_indices],
     type = "l", col = "darkgreen", lwd = 1,
     xlab = "Position (Mb)", ylab = "cM",
     main = "Genetic Map (cM) across Chromosome 1")

par(mfrow = c(1, 1))
```

## 3. Matching SNP Positions to Recombination Rates

To use recombination rates with your gene data, you need to map gene
positions to genetic coordinates.

### 3.1 Gene Midpoint Calculation

polylinkR automatically calculates gene midpoints when genomic
coordinates are provided:

``` r
# Load example gene data
data(Anatolia_EF_CLR)

# Calculate midpoints
gene_data <- copy(Anatolia_EF_CLR)
gene_data[, midpos := (startpos + endpos) / 2]

# Show first few genes with their coordinates
head(gene_data[, .(objID, objName, chr, startpos, endpos, midpos)])
```

### 3.2 Interpolating Genetic Positions

The [`plR_read()`](../reference/plR_read.md) function handles this
conversion automatically when a `rec.rate` file is provided. Here’s how
the interpolation works internally:

``` r
# Function to interpolate genetic positions
interpolate_cM <- function(chr_num, positions, rr_data) {
  rr_chr <- rr_data[chr == chr_num]
  approx(rr_chr$pos, rr_chr$cM, xout = positions, rule = 2)$y
}

# Example: Get genetic positions for first 10 genes
genes_sample <- gene_data[chr %in% 1:22][1:10]

genes_sample[, genetic_pos := interpolate_cM(chr, midpos, human_sexavg_RR_Bherer2017)]

cat("Example gene positions:\n")
print(genes_sample[, .(objName, chr, midpos, genetic_pos)])
```

### 3.3 Handling Edge Cases

- **Genes outside map boundaries**: Extrapolated using nearest available
  rate
- **Gaps in recombination map**: Linear interpolation between known
  points
- **Chromosome ends**: Rule=2 in
  [`approx()`](https://rdrr.io/r/stats/approxfun.html) extends boundary
  values

## 4. Using the `rr` Parameter in Functions

In polylinkR, recombination rate data is specified through the
`rec.rate.path` parameter in [`plR_read()`](../reference/plR_read.md).

### 4.1 Setting Up Input with Recombination Rates

``` r
# Load required data
data(PolyLinkR_SetInfo)
data(PolyLinkR_SetObj)

# Create temporary directory
temp_dir <- tempfile("polylinkR_rr_demo")
dir.create(temp_dir, recursive = TRUE)

# Write gene data (with genomic coordinates)
fwrite(Anatolia_EF_CLR, file.path(temp_dir, "ObjInfo.txt"), sep = "\t")
fwrite(PolyLinkR_SetInfo, file.path(temp_dir, "SetInfo.txt"), sep = "\t")
fwrite(PolyLinkR_SetObj, file.path(temp_dir, "SetObj.txt"), sep = "\t")

# Write recombination rate file
fwrite(human_sexavg_RR_Bherer2017, file.path(temp_dir, "RecRate.txt"), sep = "\t")

# Read with recombination rates
plr_with_rr <- plR_read(
  input.path = temp_dir,
  verbose = TRUE
)

# Check that coordinates were converted
print("Data read with recombination rate conversion:")
head(plr_with_rr$obj.info[, c("objID", "objName", "chr", "midpos")])
```

### 4.2 Available Mapping Functions

The `map.fun` parameter controls how recombination rates are converted
to genetic distances:

| Function            | Formula                     | Use Case                   |
|---------------------|-----------------------------|----------------------------|
| `Kosambi` (default) | `0.25 * log((1+2r)/(1-2r))` | Standard human genetics    |
| `Haldane`           | `-0.5 * log(1-2r)`          | No interference assumption |
| `Carter-Falconer`   | Complex                     | Mouse genetics             |
| `Morgan`            | `r`                         | 1:1 correspondence         |

Where `r` is the recombination fraction.

``` r
# Example: Compare mapping functions
recomb_fraction <- seq(0, 0.5, length.out = 100)

# Haldane mapping
haldane <- function(r) -0.5 * log(1 - 2 * r)

# Kosambi mapping
kosambi <- function(r) 0.25 * log((1 + 2 * r) / (1 - 2 * r))

# Plot comparison (excluding r = 0.5 where functions are undefined)
r_valid <- recomb_fraction[recomb_fraction < 0.49]

plot(r_valid, kosambi(r_valid), type = "l", col = "blue", lwd = 2,
     xlab = "Recombination Fraction (r)", ylab = "Genetic Distance (cM)",
     main = "Comparison of Mapping Functions")
lines(r_valid, haldane(r_valid), col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Kosambi (default)", "Haldane"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2))
```

## 5. Visualizing Linkage Signals Relative to Recombination Hotspots

Understanding how your signals relate to recombination hotspots can
provide biological insights.

### 5.1 Identifying Recombination Hotspots

``` r
# Define hotspots as regions with rate > 2 cM/Mb
rr_data <- copy(human_sexavg_RR_Bherer2017)
rr_data[, is_hotspot := rate > 2.0]

# Summarize by chromosome
hotspot_summary <- rr_data[, .(
  total_length_Mb = (max(pos) - min(pos)) / 1e6,
  avg_rate = mean(rate, na.rm = TRUE),
  max_rate = max(rate, na.rm = TRUE),
  hotspot_fraction = sum(is_hotspot, na.rm = TRUE) / .N
), by = chr]

print("Hotspot summary by chromosome:")
print(hotspot_summary[chr %in% 1:5])
```

### 5.2 Comparing Gene Density to Recombination Rate

``` r
# Analyze gene distribution relative to recombination
# Using chr 1 as example
chr1_genes <- Anatolia_EF_CLR[chr == 1]
chr1_rr <- human_sexavg_RR_Bherer2017[chr == 1]

# Bin genes by position and compare to recombination rate
bin_size <- 5e6  # 5 Mb bins
chr1_genes[, bin := floor(startpos / bin_size) * bin_size]
chr1_rr[, bin := floor(pos / bin_size) * bin_size]

gene_density <- chr1_genes[, .(n_genes = .N), by = bin]
rr_by_bin <- chr1_rr[, .(avg_rate = mean(rate, na.rm = TRUE)), by = bin]

# Merge and plot
combined <- merge(gene_density, rr_by_bin, by = "bin", all = TRUE)
combined[is.na(n_genes), n_genes := 0]
combined[is.na(avg_rate), avg_rate := 0]

par(mar = c(5, 4, 4, 4))
plot(combined$bin / 1e6, combined$n_genes, type = "h", col = "darkblue",
     xlab = "Position (Mb)", ylab = "Gene Count", lwd = 2,
     main = "Gene Density vs Recombination Rate (Chr 1)")
par(new = TRUE)
plot(combined$bin / 1e6, combined$avg_rate, type = "l", col = "red", lwd = 2,
     xlab = "", ylab = "", axes = FALSE)
axis(side = 4)
mtext("Recombination Rate (cM/Mb)", side = 4, line = 2.5, col = "red")
legend("topright", legend = c("Gene Density", "Recombination Rate"),
       col = c("darkblue", "red"), lwd = 2, lty = c(1, 1))
```

### 5.3 Significant Pathways and Recombination Landscape

When interpreting pathway results, consider:

- **Hotspot-enriched pathways**: May show stronger signals due to
  recombination-driven diversity
- **Cold region pathways**: LD may extend over longer distances,
  requiring larger `cgm.range`
- **Mixed regions**: Most realistic; signals reflect true biological
  enrichment

## 6. Differences Between Raw p-values and Recombination-Adjusted Results

### 6.1 Theoretical Expectations

| Scenario               | Without RR           | With RR                  | Expected Difference                   |
|------------------------|----------------------|--------------------------|---------------------------------------|
| Cold region pathway    | Underestimates LD    | Correct LD               | p-values become more conservative     |
| Hotspot pathway        | Overestimates LD     | Correct LD               | p-values may become less conservative |
| Long-range interaction | Assumes independence | Accounts for cM distance | May detect true long-range patterns   |

### 6.2 Practical Comparison

Here’s how to run analyses with and without recombination rates:

``` r
# Create two input directories
temp_dir_bp <- tempfile("polylinkR_bp")
temp_dir_cm <- tempfile("polylinkR_cm")
dir.create(temp_dir_bp, recursive = TRUE)
dir.create(temp_dir_cm, recursive = TRUE)

# Write common files
fwrite(Anatolia_EF_CLR, file.path(temp_dir_bp, "ObjInfo.txt"), sep = "\t")
fwrite(Anatolia_EF_CLR, file.path(temp_dir_cm, "ObjInfo.txt"), sep = "\t")
fwrite(PolyLinkR_SetInfo, file.path(temp_dir_bp, "SetInfo.txt"), sep = "\t")
fwrite(PolyLinkR_SetInfo, file.path(temp_dir_cm, "SetInfo.txt"), sep = "\t")
fwrite(PolyLinkR_SetObj, file.path(temp_dir_bp, "SetObj.txt"), sep = "\t")
fwrite(PolyLinkR_SetObj, file.path(temp_dir_cm, "SetObj.txt"), sep = "\t")

# Only cM directory gets recombination rates
fwrite(human_sexavg_RR_Bherer2017, file.path(temp_dir_cm, "RecRate.txt"), sep = "\t")

# Read both versions
plr_bp <- plR_read(input.path = temp_dir_bp, verbose = FALSE)
plr_cm <- plR_read(input.path = temp_dir_cm, verbose = FALSE)

# Compare coordinate systems
cat("Without recombination rates - midpos (bp):\n")
print(summary(plr_bp$obj.info$midpos))

cat("\nWith recombination rates - midpos (cM):\n")
print(summary(plr_cm$obj.info$midpos))
```

### 6.3 Interpreting Coordinate Differences

The conversion from bp to cM typically:

- **Reduces the dynamic range**: 250 Mb chromosome → ~200 cM
- **Compresses hotspots**: Regions with high recombination are “closer”
  in cM space
- **Expands cold regions**: Regions with low recombination are “farther”
  in cM space

``` r
# Compare distributions
par(mfrow = c(1, 2))

# Physical coordinates
hist(plr_bp$obj.info$midpos / 1e6, breaks = 50, 
     main = "Physical Distance (Mb)", xlab = "Position (Mb)",
     col = "lightblue", border = "white")

# Genetic coordinates
hist(plr_cm$obj.info$midpos, breaks = 50,
     main = "Genetic Distance (cM)", xlab = "Position (cM)",
     col = "lightgreen", border = "white")

par(mfrow = c(1, 1))
```

### 6.4 Impact on Autocorrelation Estimation

The [`plR_rescale()`](../reference/plR_rescale.md) function uses
coordinates for autocorrelation estimation. Key parameters affected:

| Parameter        | Description                    | Adjustment with RR                     |
|------------------|--------------------------------|----------------------------------------|
| `cgm.range`      | Maximum lag for autocovariance | Specify in cM (default: 2 cM ≈ 2-5 Mb) |
| `cgm.bin`        | Minimum gene pairs per bin     | May need reduction for sparse cM data  |
| Distance weights | Kernel weighting               | More accurate genetic proximity        |

### 6.5 Best Practices for Recombination Rate Usage

**When to use recombination rates:** - \[ \] You have access to
population-matched recombination maps - \[ \] Your analysis spans
diverse genomic regions (hotspots and cold regions) - \[ \] You’re
interested in fine-scale autocorrelation patterns - \[ \] Your species
has well-characterized recombination landscapes

**When physical distance may suffice:** - \[ \] Analyzing a focused
region with uniform recombination - \[ \] Computational constraints
prevent cM conversion - \[ \] Using species without available
recombination maps - \[ \] Preliminary exploratory analysis

**Quality checks:** - \[ \] Verify chromosome naming matches between
gene data and RR file - \[ \] Check for genes outside the recombination
map boundaries - \[ \] Compare cM chromosome lengths to expected genetic
map lengths - \[ \] Inspect plots for obvious conversion artifacts

## 7. Advanced Topics

### 7.1 Creating Custom Recombination Maps

If you have recombination data in a different format:

``` r
# Example: Converting from PLINK format
# PLINK format: chromosome, SNP, genetic position, physical position

convert_plink_to_hapmap <- function(plink_map) {
  # plink_map: data.table with columns CHR, SNP, cM, BP
  
  # Sort by chromosome and position
  setorder(plink_map, CHR, BP)
  
  # Calculate recombination rates
  plink_map[, rate := c(NA, diff(cM) / diff(BP) * 1e6), by = CHR]
  
  # Rename columns to HapMap format
  setnames(plink_map, 
           c("CHR", "BP", "rate", "cM"),
           c("chr", "pos", "rate", "cM"))
  
  return(plink_map[, .(chr, pos, rate, cM)])
}
```

### 7.2 Sex-Specific Recombination

The included Bherer2017 data is sex-averaged. For sex-specific analyses:

- Male recombination: More concentrated in hotspots
- Female recombination: More evenly distributed
- Consider using sex-specific maps if analyzing sex-biased traits

### 7.3 Population-Specific Maps

Recombination landscapes vary across populations:

- **European (CEU)**: Well-characterized, many hotspots shared
- **African (YRI)**: Different hotspot locations, generally higher
  recombination
- **Asian (CHB/JPT)**: Intermediate patterns

Use population-matched maps when available for your study cohort.

## Session Information

``` r
sessionInfo()

# Clean up
unlink(temp_dir, recursive = TRUE)
unlink(temp_dir_bp, recursive = TRUE)
unlink(temp_dir_cm, recursive = TRUE)
```

## References

1.  **deCODE Map**: Bherer C, et al. (2017) The recombination landscape
    in African Americans differs from that in Europeans. *Nature
    Genetics* 49(12):1760-1765.

2.  **HapMap Format**: The International HapMap Consortium (2007) A
    second generation human haplotype map of over 3.1 million SNPs.
    *Nature* 449:851-861.

3.  **Mapping Functions**: Haldane JBS (1919) The combination of linkage
    values and the calculation of distances between the loci of linked
    factors. *Journal of Genetics* 8:299-309.

4.  **Kosambi Function**: Kosambi DD (1944) The estimation of map
    distances from recombination values. *Annals of Eugenics*
    12:172-175.

5.  **Recombination Hotspots**: Myers S, et al. (2005) A fine-scale map
    of recombination rates and hotspots across the human genome.
    *Science* 310:321-324.

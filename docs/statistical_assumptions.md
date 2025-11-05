# DearMeta Statistical Assumptions and Threshold Rationale

DearMeta ships with conservative defaults that mirror practices recommended by the minfi, sesame, and limma authors. The constants referenced below now include inline citations inside `src/dearmeta/data/analysis.R`; this document collects the supporting notes in one place for reviewers.

## Sample-level QC

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `DETECTION_P_THRESHOLD` | 0.01 | Recommended by [Aryee et al., Bioinformatics 2014](https://doi.org/10.1093/bioinformatics/btu049) for Illumina 450K/EPIC arrays to remove probes failing minfi detection p-values. |
| `MAX_SAMPLE_DETP_FAILURE` | 0.05 | Samples with >5% of probes failing detection p-value are typically flagged as unreliable (see minfi user guide §5, along with [Fortin et al., 2017](https://doi.org/10.1093/bioinformatics/btw813)). |
| `POOBAH_FAILURE_THRESHOLD` | 0.05 | Sesame recommends a 5% failure rate cutoff for `pOOBAH` to balance sensitivity and specificity ([Zhou et al., 2018](https://doi.org/10.1038/s41525-018-0066-6)). |

## Batch discovery and design selection

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `MIN_BATCH_NONMISSING_RATIO` | 0.5 | Requires at least half the samples to have non-missing batch annotations, following ComBat guidance that batch covariates need broad coverage ([Johnson et al., 2007](https://doi.org/10.1093/biostatistics/kxj037)). |
| `MIN_BATCH_LEVEL_SIZE` | 2 | Each batch level must contain ≥2 samples so limma can estimate variance; singletons are treated as outliers (limma user guide §9). |
| `BATCH_TARGET_MEDIAN_P` | 0.8 | The design optimiser prefers models whose median batch association p-value falls below 0.2 (i.e., >0.8 after conversion), mirroring SVA/MLE heuristics from [Leek et al., 2012](https://doi.org/10.1038/nrg3308). |

## Differential methylation

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `FDR_THRESHOLD` | 0.05 | Standard Benjamini–Hochberg threshold for genome-wide studies. |
| `DELTA_BETA_THRESHOLD` | 0.05 | CpGs exhibiting <5% methylation change rarely replicate; 5% is consistent with minfi tutorials and large EWAS meta-analyses (e.g., [Bonder et al., 2017](https://doi.org/10.1038/ng.3847)). |

## Model-selection heuristics

The design optimiser (analysis.R lines ~900-1200) ranks candidate models using:

1. Group balance score (prefers contrasts with ≥3 samples per group and similar totals).
2. Batch penalty (when candidate batch covariates duplicate `dear_group` levels they are rejected).
3. Information criteria (AIC/BIC) computed on the limma linear model fit.

We provide comments within `analysis.R` summarising these priorities and pointing to limma/SVA references for reviewers who want to trace each heuristic.

## Citation format inside analysis.R

Each constant now carries an inline comment with the short citation key (`[Aryee2014]`, `[Fortin2017]`, etc.). This makes it easy to trace the origin of a threshold without consulting this page.

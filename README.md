# Herbaria aDNA Damage
This repository contains a set of bash and R scripts to replicate the analyses and generate the figures in:

**Porrelli, S., Fornasiero, A., Le, P.H., Yin, W., Navarrete Rodrigues, M., Mohammed, N.,
Himmelbach, A., Clarke, A.C., Stein, N., Kersey, P.J., Wing, R.A., Gutaker, R.M. (2025)
Patterns of aDNA Damage Through Time end Environments – lessons from herbarium specimens**  
d.o.i.: [XXXXXXXXXXXXXXX](https://linktodoi.com)

The workflow processes low-throughput (screening) sequencing data from 573 herbarium samples representing six plant species from two genera (*Hordeum* and *Oryza*) collected over a 220-year period (1797-2017) across the Americas and Eurasia.

## Summary
Ancient DNA degrades over time through two primary mechanisms: cytosine deamination and depurination. Understanding these degradation patterns is crucial for optimizing aDNA research and informing museum curation practices. While previous studies have focused on archaeological samples with highly variable preservation conditions, herbarium specimens offer a unique opportunity to study DNA degradation under standardized storage conditions, allowing detection of subtle temporal patterns that are often masked in archaeological contexts.

### DATA AVAILABILITY (SRA - NCBI):

Samples             | BioProject Link
------------------- | --------------------------
*Hordeum vulgare*   | [PRJNA1288534](https://example.com)
*Hordeum spontaneum*| [PRJNA1289164](https://example.com)
*Oryza rufipogon*   | [PRJNA1288425](https://example.com)
*Oryza grandiglumis*| [PRJNA1288424](https://example.com)
*Oryza latifolia*   | [PRJNA1288423](https://example.com)
*Oryza* Americas    | [PRJNA1302186](https://example.com)

## Workflow Description
This workflow consists of BASH and R scripts that process raw sequencing data and calculate damage metrics using the pipeline described in [Latorre *et. al.*, 2020](https://doi.org/10.1002/cppb.20121). Multiple approaches are used to characterize aDNA damage:

- Temporal analysis: Examines how DNA damage metrics change over time by relating fragment length, deamination patterns, and endogenous DNA content to specimen age.
- Environmental analysis: Investigates how climatic conditions at the time and location of specimen collection (temperature, precipitation, seasonality) influence DNA preservation using high-resolution climate data from CHELSA V2.1.
- Comparative analysis: Tests whether different plant genera (*Hordeum* vs *Oryza*) show distinct DNA damage patterns, potentially reflecting differences in tissue composition, growing environments, or biochemical properties.
- Variance partitioning: Quantifies the relative contributions of age, environment (temperature and precipitation), taxonomy (genus), and storage conditions (herbarium) to observed DNA damage patterns.


## Bioinformatics Pipeline (Processing Raw FASTQ Files)

- `01_Plant_aDNA_screening_prep.sh`: Sets up the environment and installs required software for running [Plant_aDNA_pipeline](https://gitlab.com/smlatorreo/plant-adna-pipeline).
- `02_Plant_aDNA_screening_main.sh`: Main pipeline to calculate aDNA damage metrics.









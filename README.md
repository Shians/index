# INtron differences to EXon (INdEX)

## Installation

```r
devtools::install_github("shians/index")
```

## Example usage

```r
library(INdEX)

exon <- readRDS(system.file("extdata/exon_dge.Rds", package = "INdEX"))
intron <- readRDS(system.file("extdata/intron_dge.Rds", package = "INdEX"))
group <- readRDS(system.file("extdata/group.Rds", package = "INdEX"))

x <- index_analysis(exon, intron, group)
plot_voom(x)
plot_lcpm_cor(x)
plot_index(x)
```

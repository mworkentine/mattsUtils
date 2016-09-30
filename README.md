# mattsUtils

A collection of helper functions for my analysis

## Microbiome helpers

### `rlog_norm_otus`  

Normalize an otu table with the `rlog` function from DESeq2

### `lfc_plot`  

Makes a plot that shows the log-fold change for OTUs

### `get_sig_taxa`  

Get significant taxa at a given p-value threshold

### `summarize_taxonomy`  

Summarize the relative amounts of taxa within a particular taxonomic rank within a higher rank

### `subset_taxa2`  

A programming safe version of `subset_taxa` from phyloseq

### `plot_OTUs`  

Make a box plot of individual OTUs

### `plot_genus_by_phylum`  

Bar plot of genera within a specific phylum

### `plot_num_seqs`  

Plot the number of sequences per sample

### `plot_bar2`  

An alternative to phyloseq's `plot_bar`

### `veganotu`  

Convert phyloseq otu table to vegan otu table

### `make_ordination_plots`  

Make ordination plots for a number of different oridination methods

### `get_richness`  

Helper function to get a nice data frame of richness estimates

### `psmelt_dplyr`  

A dplyr version of `phyloseq::psmelt`


## RNA-seq helpers

### `single_gene_plot`  

Make a bar or dot plot of single genes

### `make_heatmap`  

Make a heatmap of specific genes in a DESeq object

### `plot_diagnositics`

Make a p-value histogram and volcano plot from DESeq results

### `plot_featureCount_stats`

Plot stats output by `Rsubread::featureCounts`

### `plotPCA2`  

An alternative to `plotPCA` that provides an option to plot the 3rd component

-----

DISCLAIMER:  This package is primarily for my own use and put here for reference and reproducability purposes.  No support is provided and material is subject to change without notice - use at your own discrecion.  

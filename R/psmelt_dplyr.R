
# borrowed from http://chuckpr.github.io/blog/melt.html

psmelt_dplyr = function(physeq) {
    sd = data.frame(sample_data(physeq))
    tt = data.frame(tax_table(physeq)) %>% add_rownames("OTU")
    otu_tab = data.frame(otu_table(physeq), check.names = FALSE) %>% add_rownames("OTU")
    otu_tab %>%
        left_join(tt) %>%
        gather_("SampleID", "Abundance", setdiff(colnames(otu_tab), "OTU")) %>%
        left_join(sd)
}




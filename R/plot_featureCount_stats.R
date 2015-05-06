#' Plot stats from featureCounts output
#'
#' Takes the count summary data from a \code{\link[Rsubread]{featureCounts}}
#' object. Alternitavley you can load the counts.txt.summary file as a dataframe
#' and provide this if using the command line version of featureCounts.
#'
#' @param stats The count summary data
#' @param palette A valid RColorBrewer palette
#' @export

plot_featureCount_stats = function(stats, palette = "Set3"){
  stats %>% gather(Sample, count, -starts_with("Status"))%>%
    group_by(Sample) %>%
    mutate(total_reads = sum(count)) %>%
    group_by(Status, Sample) %>%
    mutate(percentage = count/total_reads) %>%
    ggplot(aes(x = Sample, y = percentage, fill = Status)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_brewer(palette = palette) +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      labs(y = "Percentage of total reads") +
      scale_y_continuous(labels = percent)

}

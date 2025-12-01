library(tidyverse)
library(cowplot)

setwd("data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation")
aed_file <- "results/final/assembly.all.maker.renamed.gff.AED.txt"

aed_cdf <- read_delim(aed_file, 
                      delim = "\t", 
                      col_names = TRUE,
                      show_col_types = FALSE)

colnames(aed_cdf) <- c("AED", "Cumulative_Fraction")

cat("\n", rep("=", 80), "\n")
cat("AED SCORE STATISTICS\n")
cat(rep("=", 80), "\n\n")

# Calculate the proportion in each AED bin
aed_data <- aed_cdf %>%
  arrange(AED) %>%
  mutate(
    Density = c(Cumulative_Fraction[1], diff(Cumulative_Fraction)),
    Percentage = Density * 100
  )

# Calculate summary statistics
fraction_below_0.25 <- aed_cdf$Cumulative_Fraction[aed_cdf$AED == 0.250]
fraction_below_0.5 <- aed_cdf$Cumulative_Fraction[aed_cdf$AED == 0.500]
fraction_above_0.5 <- 1 - fraction_below_0.5

cat("AED Quality Distribution:\n")
cat("Genes with AED ≤ 0.25 (excellent):", sprintf("%.2f%%\n", fraction_below_0.25 * 100))
cat("Genes with AED ≤ 0.5 (good):", sprintf("%.2f%%\n", fraction_below_0.5 * 100))
cat("Genes with AED > 0.5 (poor):", sprintf("%.2f%%\n", fraction_above_0.5 * 100))

# ============================================================================
# AED Distribution Histogram
# ============================================================================

p_aed_hist <- ggplot(aed_data, aes(x = AED, y = Percentage)) +
  geom_bar(stat = "identity", fill = "lightblue", color = "black", width = 0.023) +
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "darkgreen", linewidth = 1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "orange", linewidth = 1) +
  annotate("rect", xmin = 0, xmax = 0.25, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "darkgreen") +
  annotate("rect", xmin = 0.25, xmax = 0.5, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "orange") +
  annotate("rect", xmin = 0.5, xmax = 1, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "red") +
  annotate("text", x = 0.125, y = max(aed_data$Percentage, na.rm = TRUE) * 0.9, 
           label = "Excellent\n(≤0.25)", size = 3.5, fontface = "bold", color = "darkgreen") +
  annotate("text", x = 0.375, y = max(aed_data$Percentage, na.rm = TRUE) * 0.9, 
           label = "Good\n(0.25-0.5)", size = 3.5, fontface = "bold", color = "darkorange") +
  annotate("text", x = 0.75, y = max(aed_data$Percentage, na.rm = TRUE) * 0.9, 
           label = "Poor\n(>0.5)", size = 3.5, fontface = "bold", color = "darkred") +
  theme_cowplot() +
  labs(
    title = "Distribution of Annotation Edit Distance (AED) Scores",
    subtitle = sprintf("Good quality (AED ≤ 0.5): %.1f%% | Excellent (AED ≤ 0.25): %.1f%%",
                       fraction_below_0.5 * 100, fraction_below_0.25 * 100),
    x = "AED Score",
    y = "Percentage of Genes (%)"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.1))

ggsave(file.path(plot_dir, "07-AED_distribution_histogram.png"), 
       p_aed_hist, width = 12, height = 8, dpi = 300)
ggsave(file.path(plot_dir, "07-AED_distribution_histogram.pdf"), 
       p_aed_hist, width = 12, height = 8)

cat("\nAED histogram saved to:", file.path(plot_dir, "07-AED_distribution_histogram.pdf"), "\n")
```

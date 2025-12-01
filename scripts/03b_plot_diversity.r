library(tidyverse)
library(data.table)  # Load data.table (which has its own melt function)

data <- "/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/results/parseRM_results/assembly.fasta.mod.out.landscape.Div.Rname.tab"

rep_table <- fread(data, header = FALSE, sep = "\t")
rep_table %>% head()

colnames(rep_table) <- c("Rname", "Rclass", "Rfam", 1:50)
rep_table <- rep_table %>% filter(Rfam != "unknown")
rep_table$fam <- paste(rep_table$Rclass, rep_table$Rfam, sep = "/")

table(rep_table$fam)

# Use data.table::melt with proper syntax
rep_table.m <- melt(rep_table, 
                    id.vars = c("Rname", "Rclass", "Rfam", "fam"), 
                    measure.vars = as.character(1:50),
                    variable.name = "variable",
                    value.name = "value")

# Filter out divergence = 1
rep_table.m <- rep_table.m %>% filter(variable != "1")

# Arrange the data
rep_table.m$fam <- factor(rep_table.m$fam, levels = c(
  "LTR/Copia", "LTR/Gypsy", "DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", "DNA/Helitron",
  "MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM"
))

rep_table.m$distance <- as.numeric(rep_table.m$variable) / 100

# Remove helitrons and NA families
rep_table.m <- rep_table.m %>% 
  filter(fam != "DNA/Helitron") %>%
  filter(!is.na(fam))

# Check if there's data to plot
print(paste("Rows in filtered data:", nrow(rep_table.m)))
print(table(rep_table.m$fam))

# Create the plot
p <- ggplot(rep_table.m, aes(fill = fam, x = distance, weight = value / 1000000)) +
  geom_bar() +
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette = "Paired") +
  xlab("Distance") +
  ylab("Sequence (Mbp)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 9, hjust = 1), 
        plot.title = element_text(hjust = 0.5))


ggsave(filename = "/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/results/parseRM_results/plot_div.pdf", 
       plot = p,
       width = 10, height = 5, useDingbats = FALSE)

substitution_rate <- 8.22 * 10^-9
rep_table.m$age <- rep_table.m$distance/(2 * substitution_rate)
write.csv(rep_table.m, "/data/users/ltucker/A_thaliana_Taz-0_organization_and_annotation/results/parseRM_results/years_divereged.csv")
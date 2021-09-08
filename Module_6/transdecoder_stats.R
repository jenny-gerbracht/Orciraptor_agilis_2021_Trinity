library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)

source(file = "../config.txt")

stats_rnaspades <- read.table(paste0(mydir, "/Module_6/rnaspades_transdecoderstats.txt"),
                         row.names = 1)

stats_trinity <- read.table(paste0(mydir, "/Module_6/trinity_transdecoderstats.txt"),
                         row.names = 1)

stats <- cbind(stats_rnaspades, stats_trinity)
colnames(stats) <- c("rnaSPAdes", "Trinity")
stats <- rownames_to_column(stats, var = "type")

stats_long <- pivot_longer(stats, cols = c("rnaSPAdes", "Trinity"), names_to = "assembly", values_to = "counts")
stats_complete <- filter(stats_long, type == "Complete" | type == "Internal" | type == "5prime_partial" | type == "3prime_partial")
stats_orientation <- filter(stats_long, type == "Forward" | type == "Reverse")

stats_complete$type <- factor(stats_complete$type, 
                              levels = c("Internal", "5prime_partial", "3prime_partial", "Complete"))

pdf(file = paste0(mydir, "/Module_6/transdecoder_stats.pdf"))
ggplot(stats_complete, aes(x = assembly, y = counts, fill = type)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.4) + 
  scale_y_continuous(name = "proportion") +
  geom_text(aes(label = type), position = position_fill(vjust = 0.5), size=3.5) +
  scale_fill_brewer(palette = "Set3", guide = NULL, direction = -1) +
  theme_void()
dev.off()

ggplot(stats_orientation, aes(x = assembly, y = counts, fill = type)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.4) + 
  scale_y_continuous(name = "proportion") +
  geom_text(aes(label = type), position = position_fill(vjust = 0.5), size=3.5) +
  scale_fill_brewer(palette = "Pastel2") +
  theme_void()


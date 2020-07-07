library(tidyverse)
require(ggthemr)
library("ggpubr")

ggthemr("flat", text_size=18)

theme_set(theme_minimal())

# With FLAIR, StringTie2, Tmerge1 (fuzz=1, tolerance1, support=4), 
# tmerge2_support (fuzz=1, tolerance=1, support=4, isoform_fraction=0) and tmerge2 (fuzz=1, tolerance=1, support=1, isoform_fraction=0.055)
DATASET = "~/Documents/BIST/Major Project/data/compare/25_06_20.csv"

OUTPUT <- "~/Documents/BIST/Major Project/results/compare/"
setwd(OUTPUT)

save_plot <- function(plot, filename, width=40, height=20) {
  ggsave(filename, plot=plot, width=width, height=height, units="cm", device="png", dpi="retina")
}

# Plots from Compare pipeline
data = read_csv(
  DATASET, 
  col_names = c("source", "dataset", "time", "peak_mem", "sensitivity", "precision"),
  col_types = cols(source = col_character(), dataset = col_character(), time = col_double(), peak_mem = col_double(), sensitivity = col_double(), precision = col_double())
) %>% 
  mutate(s_and_p = (sensitivity + precision) / 2) %>%
  pivot_longer(c(sensitivity, precision), names_to="key", values_to="value") %>% 
  filter(!is.na(source) | is.na(dataset)) %>%
  # Remove malformed rows
  filter(startsWith(dataset, "ont") | startsWith(dataset, "pacBio")) %>%
  mutate(time = log(time)) %>% 
  mutate(peak_mem = log(peak_mem)) %>%
  mutate(sequencer = ifelse(str_detect(dataset, "^ont"), "ONT", "PacBio"))
  
pacbio_data <- data %>% filter(sequencer == "PacBio")
ont_data <- data %>% filter(sequencer == "ONT")

pretty_names <- scale_x_discrete(labels=c("FLAIR" = "FLAIR", "StringTie2" = "StringTie2", "tmerge1" = "tmerge1", "tmerge2" = "tmerge2", "tmerge2_splice_scoring" = "tmerge2 (Splice Scoring)"))

density_plot <- list(
  geom_violin(position=position_dodge(0.9), stat="ydensity", trim = TRUE, show.legend = TRUE),
  geom_boxplot(width = 0.15, position=position_dodge(0.9), outlier.shape=NA),
  theme(legend.title=element_blank(), text=element_text(size=20)),
  pretty_names
)

bar_plot <- list(
  geom_bar(stat="identity", position="dodge"),
  theme(legend.title=element_blank(), text=element_text(size=20)),
  pretty_names,
  scale_fill_manual(values = swatch()[3:7])
)

# Density plots for ONT
save_plot(ont_data %>% group_by(source, key) %>% filter(!is.na(value)) %>%
  ggplot(aes(x = source, y = value, fill=key)) + 
  density_plot +
  xlab("") + ylab("Value (%)"),
filename="ont_s_and_p.png")

save_plot(ont_data %>% group_by(source) %>% filter(!is.na(time)) %>%
  ggplot(aes(x = source, y = time, fill=source)) +
  density_plot +
  theme(legend.position = "none") +
  xlab("") + ylab("Time ln(s)"),
filename="ont_times.png")

save_plot(ont_data %>% group_by(source) %>% filter(!is.na(peak_mem)) %>%
  ggplot(aes(x = source, y = peak_mem, fill=source)) +
  density_plot +
  theme(legend.position = "none") +
  xlab("") + ylab("Peak memory ln(kb)"),
filename="ont_mem.png")

# Bar plots for PacBio
save_plot(pacbio_data %>% group_by(source) %>% filter(!is.na(s_and_p)) %>%
  ggplot(aes(x = source, y = s_and_p, fill=dataset)) + 
  bar_plot +
  xlab("") + ylab("[(Sn+Pn) / 2] (%)"),
filename="pb_s_and_p.png")

save_plot(pacbio_data %>% group_by(source) %>% filter(!is.na(time)) %>%
  ggplot(aes(x = source, y = time, fill=dataset)) +
  bar_plot +
  xlab("") + ylab("Time ln(s)"),
filename="pb_time.png")

save_plot(pacbio_data %>% group_by(source) %>% filter(!is.na(peak_mem)) %>%
  ggplot(aes(x = source, y = peak_mem, fill=dataset)) +
  bar_plot +
  xlab("") + ylab("Peak memory ln(kb)"),
filename="pb_mem.png")
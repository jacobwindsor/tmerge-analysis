library(tidyverse)
library("ggpubr")
library(tidyverse)
require(ggthemr)
require(interactions)
require(mgcv)
library(formattable)
library(sparkline)
ggthemr("flat", text_size=16)

#Set a few color variables to make our table more visually appealing

blueDark = "#035AA6"
blueLight = "#4192D9"
greenDarkest = "#026873"
greenDark = "#038C8C"
greenLight = "#03A696"

DATASET_SUPPORT <- ""
DATASET_ISOFORM_FRACTION <- ""

OUTPUT <- "~/Documents/BIST/Major Project/results/options_compare_support_isoform_fraction/"
setwd(OUTPUT)

save_plot <- function(plot, filename, width=20, height=20) {
  ggsave(filename, plot=plot, width=width, height=height, units="cm", device="png", dpi="retina")
}

import_data <- function(path) {
  data_raw <- read_csv(
    path, 
    col_names = c("dataset", "tolerance", "support_or_isoform_fraction", "fuzz", "sensitivity", "precision"),
    col_types = cols(dataset = col_character(), tolerance = col_double(), support_or_isoform_fraction = col_double(), fuzz = col_double(), sensitivity = col_double(), precision = col_double())
  ) %>%
    filter(str_starts(dataset, "ont") | str_starts(dataset, "pacBio")) %>% # Remove any malformed rows
    mutate(sequencer = ifelse(str_detect(dataset, "^pacBio"), "PacBio", "ONT")) %>%
    mutate(s_and_p = (sensitivity + precision) / 2)
}

# Step 1. Import data and wrangle
# ===============================
data_raw <- import_data(DATASET_ISOFORM_FRACTION) %>% 
  mutate(type = "min_isoform_fraction") %>% 
  full_join(import_data(DATASET_SUPPORT) %>% mutate(type = "min_read_support"))

data_raw <- data_raw %>% 
  drop_na() %>%
  group_by(dataset, type)

data <- data_raw %>% filter(sequencer == "ONT") # Only do for ONT as not enough datasets from PacBio

# STEP 2. Compare variance between min_read_support and min_isoform_fraction optimum values
# =========================================================================================
save_plot(
  data %>%
    filter(tolerance==0, fuzz==0) %>% 
    filter(sensitivity > 25, s_and_p == max(s_and_p)) %>% 
    filter(support_or_isoform_fraction == min(support_or_isoform_fraction)) %>% 
    ggplot(aes(x=type, y=support_or_isoform_fraction, fill=type)) + 
    ylab("Optimum Value") +
    geom_boxplot(width = 0.15, position=position_dodge(0.9), outlier.shape=NA, show.legend = FALSE),
  filename="variance_isoform_fraction_vs_support.png")

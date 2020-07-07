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

# Description:
# From options pipeline, find the most optimium value for each option

DATASET <- "~/Documents/BIST/Major Project/data/options/22_06_20.csv"

OUTPUT <- "~/Documents/BIST/Major Project/results/options_support/"
setwd(OUTPUT)

save_plot <- function(plot, filename, width=20, height=20) {
  ggsave(filename, plot=plot, width=width, height=height, units="cm", device="png", dpi="retina")
}

import_data <- function(path) {
  data_raw <- read_csv(
    path, 
    col_names = c("dataset", "tolerance", "support", "fuzz", "sensitivity", "precision"),
    col_types = cols(dataset = col_character(), tolerance = col_double(), support = col_double(), fuzz = col_double(), sensitivity = col_double(), precision = col_double())
  ) %>%
    filter(str_starts(dataset, "ont") | str_starts(dataset, "pacBio")) %>% # Remove any malformed rows
    mutate(sequencer = ifelse(str_detect(dataset, "^pacBio"), "PacBio", "ONT")) %>%
    mutate(s_and_p = (sensitivity + precision) / 2)
}

# Step 1. Import data and wrangle``
# ===============================

# Find the most optimum values for each option
data_raw <- import_data(DATASET)

# Print the counts for each
formattable(tribble(
  ~Type, ~N,
  "ONT",  data_raw %>% filter(sequencer == "ONT") %>% distinct(dataset) %>% count() %>% pull(n),
  "ONT SIRVs", data_raw %>% filter(sequencer == "ONT") %>% drop_na() %>% distinct(dataset) %>% count() %>% pull(n),
  "PacBio", data_raw %>% filter(sequencer == "PacBio") %>% distinct(dataset) %>% count() %>% pull(n),
  "PacBio SIRVs", data_raw %>% filter(sequencer == "PacBio") %>% drop_na() %>% distinct(dataset) %>% count() %>% pull(n)
))

data_raw <- data_raw %>% 
  drop_na() %>%
  group_by(dataset)

data <- data_raw %>% filter(sequencer == "ONT") # Only do for ONT as not enough datasets from PacBio

plot_option <- function(data, option_name) {
  data %>%
    pivot_longer(c(tolerance, support, fuzz), names_to="option", values_to="option_value") %>%
    filter(option==option_name) %>%
    ungroup() %>% group_by(option_value) %>%
    pivot_longer(c(sensitivity, precision), names_to="y_type", values_to="y_value") %>%
    ggplot(aes(x=option_value, y=y_value, linetype=y_type, color=y_type, fill=y_type)) +
    ylab("%") + xlab("Option Value") + ylim(0, 100) +
    theme(legend.title=element_blank(), plot.margin = margin(1,0.2,1,0.2, "cm")) +
    geom_smooth()
}

# Step 2. Plot each option changing while holding others
# TODO: get these into one graph
# ======================================================

p_fuzz <- data %>%
  filter(support==4,tolerance==0) %>%
  plot_option(option_name="fuzz") + xlab("end_fuzz (nt)")

p_support <- data %>%
  filter(fuzz==0,tolerance==0) %>%
  plot_option(option_name="support") +
  scale_x_continuous(breaks = c(2,4,6,8,10)) + xlab("min_read_support")
  
p_tolerance <- data %>%
  filter(support==1,fuzz==0) %>%
  plot_option(option_name="tolerance") + xlab("exon_overhang_tolerance (nt)")

save_plot(p_fuzz, filename="fuzz.png")
save_plot(p_support, filename="support.png")
save_plot(p_tolerance, filename="tolerance.png")

ggexport(
  ggarrange(
    p_fuzz, p_support, p_tolerance,
    labels=c("end_fuzz", "min_read_support", "tolerance"),
    font.label = list(size = 11, color = "black", face = "bold", family = NULL),
    vjust = 2,
    common.legend = TRUE, legend = "right",
    nrow=2, ncol=2
  ), 
  filename="options.png" 
)

# Step 3: Find interactions
# From https://m-clark.github.io/generalized-additive-models/application.html#multiple-predictors
# and http://environmentalcomputing.net/intro-to-gams/
# =========================

# make a grid of possible tolerances and supports
poss_tolerances <- seq(0,20)
poss_supports <- seq(1,10)
poss_fuzzs <- seq(0,5)
grid_support <- expand.grid(tolerance = poss_tolerances, fuzz = poss_fuzzs, support = poss_supports)



# Model support
# ================
# Tolerance affects precision but not sensitivity
mod_sensitivity = gam(sensitivity ~ s(support, k = 6) + s(fuzz, k=5), data=data)
mod_precision = gam(precision ~ s(tolerance) + s(support, k = 6) + s(fuzz, k=5), data=data)

predict_sensitivity <- predict(mod_sensitivity, grid_support)
predict_precision <- predict(mod_precision, grid_support)

predict_df <- grid_support %>%
  mutate(sensitivity = predict_sensitivity) %>%
  mutate(precision = predict_precision) %>%
  mutate(s_and_p = (sensitivity + precision) / 2)

# Make a histogram of sensitivity and precisions to decide cutoff
save_plot(
  predict_df %>% 
    pivot_longer(c(sensitivity, precision), names_to="key", values_to="value") %>% 
    ggplot(aes(x=value, fill=key)) + 
    geom_histogram(alpha=0.5, position="identity") +
    xlab("%") +
    scale_x_continuous(breaks = seq(0, 100, by = 10)) +
    theme(legend.title=element_blank()),
  filename = "ont_predicted_support_hist.png",
  width = 30
)

pred_df <- predict_df %>%
  filter(sensitivity > 40 & sensitivity > 40)

print(pred_df %>% filter(s_and_p == max(s_and_p)))
formattable(pred_df %>% filter(s_and_p == max(s_and_p)))

save_plot(
  predict_df %>% filter(fuzz == 2) %>%
    ggplot() +
    aes(x = tolerance, y = support, z = s_and_p) + 
    xlab("tolerance") + ylab("min_read_support") +
    geom_contour_filled(),
  filename="contour_support_fuzz-0.png")

save_plot(
  predict_df %>% filter(tolerance == 0) %>%
    ggplot() +
    aes(x = fuzz, y = support, z = s_and_p) + 
    xlab("end_fuzz") + ylab("min_read_support") +
    geom_contour_filled(),
  filename="contour_support_tolerance-0.png")

# Step 5. analysis for PacBio
# ===========================
data <- data_raw %>% filter(sequencer == "PacBio") %>% drop_na() %>% 
  group_by(dataset)

# Make a histogram of sensitivity and precisions to decide cutoff for sensitivity
save_plot(
  data %>% 
    pivot_longer(c(sensitivity, precision), names_to="key", values_to="value") %>% 
    ggplot(aes(x=value, fill=key)) + 
    geom_histogram(alpha=0.5, position="identity") +
    scale_x_continuous(breaks = seq(0, 100, by = 10)) +
    theme(legend.title=element_blank()),
  filename = "pacbio_hist.png"
)

make_summary <- function(data) {
  data %>% select(tolerance, support, fuzz) %>% 
    group_map(
      ~ .x %>% map(~ if(length(unique(.x)) > 1) { stop("Can only make summary if one option per group") })
    )
  
  data %>% summarise(
    exon_overhang_tolerance = min(tolerance),
    fuzz = min(fuzz),
    min_read_support = min(support),
    precision = max(precision),
    sensitivity = max(sensitivity)
  )
}

make_table <- function(data) {
  data %>% 
    formattable(list(
      area(col = 7:8) ~  color_tile(greenLight, greenDark),
      `precision` = color_tile(blueLight, blueDark)
    ))
}

filter_highest_s_and_p <- function(data) { data %>% filter(sensitivity > 60) %>% filter(s_and_p == max(s_and_p)) %>% filter(tolerance == min(tolerance))}

# Filter out precison > 80% and sensitivyt < 20% as shown in histogram
data %>% filter_highest_s_and_p() %>% mutate(method = "min_read_support") %>% make_summary() %>% make_table()

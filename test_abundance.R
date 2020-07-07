# Find the most optimum values for each option
data = read_csv(
  "~/Documents/BIST/Major Project/data/abundance/b.csv", 
  col_names = c("dataset", "tolerance", "support", "abundance", "fuzz", "sensitivity", "precision"),
  col_types = cols(dataset = col_character(), tolerance = col_double(), support = col_double(), abundance = col_double(), fuzz = col_double(), sensitivity = col_double(), precision = col_double())
) %>%
  drop_na()


data %>% 
  mutate(s_and_p = sensitivity + precision) %>%
  group_by(dataset) %>%
  filter(s_and_p == max(s_and_p)) %>%
  summary(mean = mean(value))

# Plot a distribution for the best options
data %>% 
  filter(abundance==0.15) %>% 
  pivot_longer(c(sensitivity, precision), names_to="key", values_to="value") %>%
  ggplot(aes(x = abundance, y = value, fill=key)) + geom_violin(position=position_dodge(0.9), stat="ydensity", trim = TRUE, show.legend = legend) +
  geom_boxplot(width = 0.15, position=position_dodge(0.9), outlier.shape=NA) +
  ggtitle("Precision and sensitivity values between tools") +
  xlab("") + ylab("Value (%)")
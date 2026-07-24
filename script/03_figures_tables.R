# =============================================================================
# 03_figures_tables.R
# -----------------------------------------------------------------------------
# Assemblage-comparability outputs for Section 3.1:
#   Table 2  — blade/bladelet counts per assemblage (includes RB_A1)
#   Figure 2 — width distributions of complete laminar blanks
#   Figure 3 — length vs width of blanks (by cortex) vs last core removals
#
# Reads the RAW attribute data directly (independent of 01/02); no upstream
# .rds needed. Uses a different library stack from the statistical scripts.
#
# INPUTS : data/dataset_blanks.csv, data/dataset_cores.csv,
#          data/dataset-for-tables.csv
# WRITES : output/figures/Figure_2.tif, Figure_3.tif  (Table 2 printed inline)
# =============================================================================

library(tidyverse)
library(janitor)
library(flextable)
library(forcats)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# Shared display order (RB_A1 kept here; it is excluded from the FLEXDIST runs).
layer_order      <- c("Pa_24A", "Ca_Au", "CTC_gic-ars", "CTC_rsa'",
                      "RB_A2", "RB_A1", "RF_A2-A1")
layer_order_disp <- gsub("_", ", ", layer_order)   # "Pa, 24A" etc. for facets

# ---- Data -------------------------------------------------------------------
dataset_blanks <- read.csv("data/dataset_blanks.csv") %>%
  filter(rowSums(is.na(across(Length:Distal.end))) <= 2)

dataset_cores <- read.csv("data/dataset_cores.csv") %>%
  filter(Layer != "RF_D3b.alpha")

dataset_blanks_tables <- read.csv("data/dataset-for-tables.csv") %>%
  filter(rowSums(is.na(across(Length:Distal.end))) <= 2)


# =============================================================================
# TABLE 2 — blank-type counts per assemblage
# =============================================================================
Table_2 <- dataset_blanks_tables %>%
  filter(Layer != "RB_A1") %>%
  mutate(Layer = fct_relevel(Layer,
           "Pa_24A", "Ca_Au", "CTC_gic-ars", "CTC_rsa'", "RB_A2", "RF_A2-A1")) %>%
  tabyl(Blank.type, Layer) %>%
  adorn_totals("col") %>% adorn_totals("row") %>%
  adorn_percentages("col") %>% adorn_pct_formatting(digits = 1) %>%
  adorn_ns(position = "front") %>%
  flextable() %>% theme_vanilla() %>%
  set_table_properties(layout = "autofit") %>% autofit()

Table_2


# =============================================================================
# FIGURE 2 — width distributions of complete laminar blanks
# =============================================================================
filtered_data <- dataset_blanks %>%
  filter(Blank.type %in% c("Blade", "Bladelet"),
         Preservation == "Complete",
         Class == "Blank")

# Balance to the smallest assemblage so distributions are visually comparable.
set.seed(123)
min_size <- filtered_data %>% count(Layer) %>% summarise(m = min(n)) %>% pull(m)

balanced_data <- filtered_data %>%
  group_by(Layer) %>% slice_sample(n = min_size) %>% ungroup() %>%
  mutate(Layer = gsub("_", ", ", Layer),
         Layer = factor(Layer, levels = layer_order_disp))

Figure_2 <- ggplot(balanced_data, aes(x = Width)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0,
                 fill = "#88CCEE", colour = "grey20", linewidth = 0.25) +
  geom_density(colour = "#332288", linewidth = 1.1) +
  geom_vline(xintercept = 12, linetype = "dashed", colour = "#CC3311",
             linewidth = 0.7) +
  annotate("text", x = 12, y = Inf, label = "12 mm",
           hjust = -0.1, vjust = 1.4, size = 3, colour = "#CC3311") +
  facet_wrap(~ Layer, scales = "free_y", ncol = 4) +
  labs(x = "Width (mm)", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey85", colour = NA),
        strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

Figure_2
ggsave("output/figures/Figure_2.tif", Figure_2, width = 10, height = 7,
       dpi = 300, units = "in", compression = "lzw")


# =============================================================================
# FIGURE 3 — blank length vs width (by cortex) against last core removals
# =============================================================================
optimal_blanks <- filtered_data %>%
  filter(Cortex %in% c("0%", "1-33%"), Preservation == "Complete") %>%
  select(Length, Width, Layer)

cortical_blanks <- filtered_data %>%
  filter(Cortex %in% c("33-66%", "66-99%", "100%"), Preservation == "Complete") %>%
  select(Length, Width, Layer)

last_removals <- dataset_cores %>%
  filter(Last.removal.blank %in% c("Blade", "Bladelet")) %>%
  select(Length = Last.removal.length, Width = Last.removal.width, Layer)

cat_levels <- c("Blank <33% cortex", "Blank >33% cortex", "Last laminar removal")

fig3_data <- bind_rows(
  mutate(optimal_blanks,  Category = cat_levels[1]),
  mutate(cortical_blanks, Category = cat_levels[2]),
  mutate(last_removals,   Category = cat_levels[3])) %>%
  mutate(Layer    = gsub("_", ", ", Layer),
         Layer    = factor(Layer, levels = layer_order_disp),
         Category = factor(Category, levels = cat_levels))

Figure_3 <- ggplot(fig3_data,
                   aes(Length, Width, colour = Category,
                       shape = Category, linetype = Category)) +
  geom_point(alpha = 0.5, size = 1.5) +
  stat_ellipse(level = 0.95, linewidth = 0.8) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  labs(x = "Length (mm)", y = "Width (mm)",
       colour = "Category", shape = "Category", linetype = "Category") +
  facet_wrap(~ Layer) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(colour = "grey70", fill = NA),
        panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

Figure_3
ggsave("output/figures/Figure_3.tif", Figure_3, width = 10, height = 7,
       dpi = 600, units = "in", compression = "lzw")

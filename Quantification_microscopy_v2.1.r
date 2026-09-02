# Put the path to the file you want to analyze here:
path <- '/Users/Timo/Downloads/test_r-script/Results.csv'

# set threshold for size, otherwise keep them 0 or inf (unit is pixels^2)
min_area <- 400
max_area <- Inf
min_foci_size <- 3
remove_first_characters <- 28
remove_last_characters <- -13
plot_width <- 32
plot_height <- 24

conditions <- c(
  "WT-HS",
  "WT+HS",
  "WT+HS_t60",
  "WT+HS_t120",
  "F436A-HS",
  "F436A+HS",
  "F436A+HS_t60",
  "F436A+HS_t120"
)


# set threshold for relative foci intensity between 0 and 1 (1 being 100%)
# if you don't want them filtered, keep the thresholds at 0 and 1 respectively
min_intensity <- 0
max_intensity <- 1



########################################

# load required packages
if(!require(tidyverse)){
  install.packages("tidyverse")
}
if(!require(writexl)){
  install.packages("writexl")
}
if(!require(stringr)){
  install.packages("stringr")
}
library(stringr)
library(tidyverse)
library(writexl)
library(dplyr)
library(ggplot2)

# ----------------------------- Import data -----------------------------
df <- read_csv(path)

filtered_df <- df %>%
  select(Image, Cell_number, Cell_Area, Cell_Intensity, Area, IntDen) %>%
  filter(Cell_Area >= min_area & Cell_Area <= max_area) %>%
  rename(Foci_Area = Area, Foci_Intensity = IntDen) %>%
  mutate(
    Foci_Area = if_else(Foci_Area < min_foci_size, NA_real_, Foci_Area),
    Foci_Intensity = if_else(Foci_Area < min_foci_size, NA_real_, Foci_Intensity)
  )

colnames(filtered_df) <- tolower(colnames(filtered_df))
# -------------------------------------------------------------------------

# ----------------------------- Summarize per cell -----------------------------
cell_summary <- filtered_df %>%
  group_by(image, cell_number, cell_area, cell_intensity) %>%
  summarise(
    n_foci = sum(!is.na(foci_area) & foci_area > 0),
    total_foci_area = sum(foci_area, na.rm = TRUE),
    total_foci_intensity = sum(foci_intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_foci = replace_na(n_foci, 0),
    total_foci_area = replace_na(total_foci_area, 0),
    total_foci_intensity = replace_na(total_foci_intensity, 0)
  )

cell_summary <- cell_summary %>%
  mutate(relative_foci_intensity = (total_foci_intensity / cell_intensity) * 100)
# -------------------------------------------------------------------------

# ----------------------------- Prepare Prism tables -----------------------------
condition_pattern <- conditions |>
  sort(decreasing = TRUE) |>
  stringr::str_replace_all("\\+", "\\\\+") |>
  paste(collapse = "|")


cell_summary_name_shortened <- cell_summary %>%
  mutate(
    condition = stringr::str_extract(image, condition_pattern)
  ) %>%
  select(condition, n_foci, relative_foci_intensity)


if (any(is.na(cell_summary_name_shortened$condition))) {
  missing_images <- cell_summary %>%
    mutate(condition = stringr::str_extract(image, condition_pattern)) %>%
    filter(is.na(condition)) %>%
    distinct(image)
  
  stop(
    "Some images could not be matched to a condition:\n",
    paste(missing_images$image, collapse = "\n")
  )
}


# Table 1: number of foci per cell
foci_count_table <- cell_summary_name_shortened %>%
  select(condition, n_foci) %>%
  rename(Foci_count = n_foci) %>%
  group_by(condition) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = condition,
    values_from = Foci_count
  ) %>%
  select(-row_id)

# Table 2: relative foci intensity per cell
relative_intensity_table <- cell_summary_name_shortened %>%
  select(condition, relative_foci_intensity) %>%
  group_by(condition) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = condition,
    values_from = relative_foci_intensity
  ) %>%
  select(-row_id)
# -------------------------------------------------------------------------

# ----------------------------- Create binned foci table -----------------------------
foci_count_table_binned <- cell_summary_name_shortened %>%
  mutate(Foci_Bin = if_else(n_foci >= 9, 9L, n_foci))

# Summarize binned counts per condition
foci_count_binned_summary <- foci_count_table_binned %>%
  group_by(condition, Foci_Bin) %>%
  summarise(n_cells = n(), .groups = "drop")
# -------------------------------------------------------------------------

# ----------------------------- Export Excel -----------------------------
new_path <- str_replace(path, ".csv", "_summary.xlsx")
output_tables <- list(
  "Filtered_Data" = cell_summary,
  "Foci_Count_Summary" = foci_count_table,
  "Relative_Intensity_Summary" = relative_intensity_table,
  "Foci_Count_Binned" = foci_count_binned_summary
)

write_xlsx(output_tables, new_path)
# -------------------------------------------------------------------------

# ----------------------------- Plotting -----------------------------
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
bwPalette <- c("#ffffff", "#b9b9b9", "#8e8e8e", "#727272", "#000000", "#0072B2", "#D55E00", "#CC79A7")

my_favourite_theme <- theme(
  panel.grid = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 14),
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 12),
  axis.text.x = element_text(angle = 45, hjust = 1)
)

folder_path <- dirname(path)

# ---- Foci count stacked bar plot ----
max_foci <- max(foci_count_table_binned$Foci_Bin, na.rm = TRUE)
foci_levels <- 0:max_foci
foci_labels <- ifelse(foci_levels == 9, "≥9", as.character(foci_levels))
gradient_colors <- colorRampPalette(c("white", "black"))(length(foci_levels))

stacked_df <- foci_count_table_binned %>%
  group_by(condition, Foci_Bin) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(
    Foci_Bin = factor(Foci_Bin, levels = foci_levels, labels = foci_labels)
  )

# ---- Foci count stacked bar plot ----
p_foci_bar <- ggplot(stacked_df, aes(x = condition, y = prop, fill = Foci_Bin)) +
  geom_bar(stat = "identity", color = "black", position = position_stack(reverse = TRUE), width = 0.5) +
  my_favourite_theme +
  scale_fill_manual(values = gradient_colors, name = "Foci per cell", guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0,0))) +
  labs(x = "Condition", y = "Proportion of cells", title = "Normalized distribution of foci counts per condition")

ggsave(file.path(folder_path, "foci-count_barplot.pdf"), plot = p_foci_bar, width = plot_width, height = plot_height, units = "cm")

# ---- Foci count violin plot with points ----
y_breaks <- seq(0, max_foci + 1, by = 2)

p_foci_violin <- ggplot(cell_summary_name_shortened, aes(x = condition, y = n_foci)) +
  geom_violin(fill = "orange", color = "black", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.6, color = "black") +
  my_favourite_theme +
  scale_y_continuous(
    breaks = y_breaks,  # discrete steps of 2
    minor_breaks = NULL                 # optional: remove minor breaks
  ) +
  labs(x = "Condition", y = "Number of Foci per Cell", title = "Distribution of Number of Foci per Cell")

ggsave(file.path(folder_path, "foci-count_violin_with_points.pdf"), plot = p_foci_violin, width = plot_width, height = plot_height, units = "cm")

# ---- Relative intensity plot with SD and points ----

# ---- Relative intensity statistics ----
relative_intensity_stats <- cell_summary_name_shortened %>%
  group_by(condition) %>%
  summarise(
    mean_relative_intensity = mean(relative_foci_intensity, na.rm = TRUE),
    sd_relative_intensity   = sd(relative_foci_intensity, na.rm = TRUE),
    n_cells                 = sum(!is.na(relative_foci_intensity)),
    .groups = "drop"
  )

# ---- Compute shared label positions ----
label_positions <- cell_summary_name_shortened %>%
  group_by(condition) %>%
  summarise(
    max_point = max(relative_foci_intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(relative_intensity_stats, by = "condition") %>%
  mutate(
    label_y = pmax(
      max_point,
      mean_relative_intensity + sd_relative_intensity,
      na.rm = TRUE
    ) * 1.05
  )


p_intensity <- ggplot(
  relative_intensity_stats,
  aes(x = condition, y = mean_relative_intensity)
) +
  # SD error bars (thick)
  geom_errorbar(
    aes(
      ymin = mean_relative_intensity - sd_relative_intensity,
      ymax = mean_relative_intensity + sd_relative_intensity
    ),
    width = 0.25,
    linewidth = 0.8
  ) +
  # n labels
  geom_text(
    data = label_positions,
    aes(
      x = condition,
      y = label_y,
      label = paste0("n=", n_cells)
    ),
    size = 3
  ) +
  # Individual data points
  geom_jitter(
    data = cell_summary_name_shortened,
    aes(x = condition, y = relative_foci_intensity),
    width = 0.15,
    height = 0,
    size = 1.2,
    alpha = 0.6,
    color = "black"
  ) +
  # Mean as thick horizontal line
  geom_point(
    size = 24,
    shape = 95,   # horizontal line
    stroke = 3,
    color = "orange"
  ) +
  my_favourite_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Condition",
    y = "Relative Foci Intensity (%)",
    title = "Mean ± SD Relative Foci Intensity per Condition"
  )


ggsave(file.path(folder_path, "relative_foci_intensity_barplot_with_SD.pdf"), plot = p_intensity, width = plot_width, height = plot_height, units = "cm")

# ---- Relative intensity violin plot with points ----

p_violin <- ggplot(
  cell_summary_name_shortened,
  aes(x = condition, y = relative_foci_intensity)
) +
  # geom_violin(fill = "steelblue", color = "black", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_text(
    data = label_positions,
    aes(
      x = condition,
      y = label_y,
      label = paste0("n=", n_cells)
    ),
    size = 3
  )+
  geom_jitter(
    width = 0.15,
    height = 0,
    size = 1,
    alpha = 0.6,
    color = "black"
  ) +
  my_favourite_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Condition",
    y = "Relative Foci Intensity (%)",
    title = "Distribution of Relative Foci Intensity per Condition"
  )


ggsave(file.path(folder_path, "relative_foci_intensity_violin_with_points.pdf"), plot = p_violin, width = plot_width, height = plot_height, units = "cm")


# ----------------------------- Multi-panel figure -----------------------------
if(!require(patchwork)) install.packages("patchwork")
library(patchwork)

# Arrange the four plots in 2x2 grid
multi_panel <- (p_foci_bar | p_foci_violin) / (p_intensity | p_violin)

# Optional: add labels to panels (A, B, C, D)
multi_panel <- multi_panel + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Save multi-panel figure as PDF
multi_panel_filename <- file.path(folder_path, "multi_panel_foci_analysis.pdf")

ggsave(
  filename = multi_panel_filename,
  plot = multi_panel,
  width = plot_width * 2,   # wider to accommodate 2 columns
  height = plot_height * 2, # taller to accommodate 2 rows
  units = "cm"
)

library(tidyverse)
library(tidyplots)

# TODO The following directory must contain multiple CSV files.
# TODO Each CSV file must contain the columns: sample, variant_count, Source.
# TODO The sample column should contain sample IDs, variant_count is the variants per genome for each individual, and Source contains the group/category.
variants_id_dir <- '../output/subset/'
# Read all CSV files in the variant directory
# Merge them into a single data frame
variants_id_files <- list.files(variants_id_dir, pattern = '\\.csv$', full.names = TRUE)

variants_id_files <- variants_id_files[file.exists(variants_id_files)]
variants_id_list <- lapply(variants_id_files, function(file) {
  read.csv(file, header = TRUE, stringsAsFactors = FALSE)
})

variants_id_df <- do.call(rbind, variants_id_list)
variants_id_df |> as_tibble() -> df_variants_id

df_variants_id |>
   mutate(variant_count = variant_count / 1000) -> df_variants_id
df_variants_id |> select(Source) |> distinct() |> pull() 

df_variants_id |> 
    tidyplot(x = Source, y = variant_count, color = Source) |>
    add_data_points_jitter(
    alpha = 0.2,
    size = 0.1,
    jitter_width = 5,
    jitter_height = 7,
    dodge_width = 1,
    cex = 1,
    preserve = "single",
    rasterize = TRUE,
    rasterize_dpi = 300
    ) |>
    add_boxplot(
      alpha = 0.8,
      show_outliers = FALSE) |>
    add_median_value(hjust = 0.5, vjust = -3) |>  
    adjust_y_axis_title('Variants per genome (K)') |>
    adjust_x_axis_title('Continent') |>
    adjust_x_axis(rotate_labels = 90) |>
    sort_x_axis_labels(.reverse = TRUE) |>
    save_plot('Boxplot.pdf')

df_variants_id |> 
    tidyplot(x = Source, y = variant_count, color = Source) |>
    add_data_points_beeswarm(
    alpha = 0.1,
    size = 0.1,
    jitter_width = 5,
    jitter_height = 7,
    dodge_width = 1,
    corral = "wrap",
    corral.width = 0.5,
    cex = 1,
    preserve = "single",
    rasterize = TRUE,
    rasterize_dpi = 600
    ) |>
    add_boxplot(
      alpha = 0.8,
      show_outliers = FALSE) |>
    add_median_value(hjust = 0.5, vjust = -3) |>  
    adjust_y_axis_title('Variants per genome (K)') |>
    adjust_x_axis_title('Continent') |>
    adjust_x_axis(rotate_labels = 90) |>
    sort_x_axis_labels(.reverse = TRUE) |>
    save_plot('Boxplot.pdf')


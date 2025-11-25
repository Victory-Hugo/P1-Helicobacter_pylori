library(tidyverse)
library(tidyplots)

df_variants <- read.csv('../output/compare_eas_global/csv/Global.csv', sep = ',')

# Convert to tibble
df_variants |> as_tibble() -> df_variants
df_variants
#* First, draw the overall bar plot
df_variants|> 
   filter(Category != 'Type') |>
   filter(Class != '') |> # Filter out rows with empty Class
   mutate(Cat = as.factor(Category)) |>
   tidyplot(x = Class, y = Count, fill = Class) |>
   add_sum_bar() |> # Add total bar
   add_sum_value(fontsize = 4) |> # Add value labels
   adjust_x_axis(rotate_labels = 90) |> # Adjust x-axis labels
   sort_x_axis_labels(by = Cat) |> # Sort by Category
   save_plot('overall_variant_statistics.pdf') 

#* Then, draw grouped bar plots
df_variants |> 
   filter(Category != 'Type') |> 
   filter(Class != '') |> # Filter out rows with empty Class
   mutate(Cat = as.factor(Category)) |>
   tidyplot(x = Class, y = Count, fill = Class) |>
   add_sum_bar() |> # Add total bar
   add_sum_value(fontsize = 4) |> # Add value labels
   adjust_x_axis(rotate_labels = 90) |> # Adjust x-axis labels
   sort_x_axis_labels(by = Cat) |> # Sort by Cat
   split_plot(by = Source, ncol = 2,widths =50,heights = 30) |>
   save_plot('variant_statistics_by_group.pdf') 

#* Draw bar plot for SNP and INDEL
df_variants |> 
   filter(Category == 'Type') |>
   filter(Class != '') |> # Filter out rows with empty Class
   mutate(Cat = as.factor(Category)) |>
   tidyplot(x = Class, y = Count, fill = Class) |>
   add_sum_bar() |> # Add total bar
   add_sum_value(fontsize = 4) |> # Add value labels
   sort_x_axis_labels(by = Cat) |> # Sort by Category
   save_plot('overall_variant_statistics_SNP_INDEL.pdf') 

#* Then, draw grouped bar plots for SNP and INDEL
df_variants |> 
   filter(Category == 'Type') |> 
   filter(Class != '') |> # Filter out rows with empty Class
   mutate(Cat = as.factor(Category)) |>
   tidyplot(x = Class, y = Count, fill = Class) |>
   add_sum_bar() |> # Add total bar
   add_sum_value(fontsize = 4) |> # Add value labels
   adjust_x_axis(rotate_labels = 90) |> # Adjust x-axis labels
   sort_x_axis_labels(by = Cat) |> # Sort by Category
   split_plot(by = Source, ncol = 1,widths =50,heights = 30) |>
   save_plot('variant_statistics_SNP_INDEL_by_group.pdf') 

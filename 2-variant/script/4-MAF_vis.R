library(tidyverse)
library(tidyplots)

df_MAF <- read.csv('../../1-Variants-stat/output/Bin_MAF_Comparison.csv', sep = ',')
df_MAF |> as_tibble() -> df_MAF

df_MAF |> 
  pivot_longer(cols = -MAF, names_to = 'Source', values_to = 'global_count') |>
  filter(Source != 'global_count') -> df_long 

df_long |> select(MAF) |> distinct() |> pull(MAF) -> MAF
MAF

maf_order <- c('<0.1', '0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-1.0',
               '1-2', '2-5', '5-10', '10-20', '20-30', '30-40', '40-50')

df_long |> 
  mutate(MAF = factor(MAF, levels = maf_order)) -> df_long
print(df_long,n = 50)

df_long  |> 
  tidyplot(x = MAF, y = global_count, fill = Source) |>
  add_barstack_absolute()  |>
  add_data_labels(label = global_count, position = "stack", size = 3) |>
  adjust_x_axis(rotate_labels = 90) |>
  save_plot("MAF_distribution.pdf")

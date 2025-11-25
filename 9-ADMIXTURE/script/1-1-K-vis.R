library(tidyverse)

# 读取CSV文件
K_CV <- read.csv('cv_results.csv')
K_CV
# 方法1: 保持K为因子
K_CV |> 
  ggplot(aes(x = K, y = CV, group = 1)) +  # 添加group=1使所有点连接在一起
  geom_line() +
  geom_point(size = 3) +  # 增大点的大小以便更清晰
  labs(title = 'Cross-validation Results for K', 
       x = 'K', 
       y = 'Cross-Validation Error') +
  scale_x_continuous(breaks = K_CV$K) +  # 确保X轴只显示数据中的K值
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) 

library(tidyplots)
K_CV |> 
   tidyplot(
    x = K,
    y = CV,
    color= CV
   ) |>
   add_line(group = 1) |>
   add_data_points(size = 3) |>
   add_data_labels(label = CV,fontsize = 7) |>
   adjust_x_axis(breaks = K_CV$K)  |>
   adjust_colors(colors_continuous_cividis) |>
   save_plot("最优K可视化.pdf")

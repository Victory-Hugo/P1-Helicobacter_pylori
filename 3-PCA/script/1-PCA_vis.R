# ============================================================================
# PCA scatter plot visualization
# Using base R plotting functions with a classic PCA style
# ============================================================================

# Load required packages
library(tidyverse)
library(RColorBrewer)
# ============================================================================
# 
# ============================================================================
PCA_result_file <- "PCA.csv" #! Must contain columns: ID, PC1, PC2, Class_big, Class_small
color_file <- "color.csv" #! Must contain columns: Class_big and color
OUT_pdf_file <- "PCA.pdf" 


# Read PCA result data
cat("Reading PCA data...\n")
# !Note: modify the path according to your actual data location
data <- read.csv(PCA_result_file, header = TRUE)

# Basic data information
cat("Data overview:\n")
cat("Number of samples:", nrow(data), "\n")
cat("Number of variables:", ncol(data), "\n")
print(str(data))

# Extract classification information
Class_big_levels <- unique(data$Class_big) #todo Make sure the column names Class_big and Class_small match the CSV file
class_small_levels <- unique(data$Class_small) #todo Make sure the column names Class_big and Class_small match the CSV file

cat("\nClass_big categories (", length(Class_big_levels), "):", paste(Class_big_levels, collapse=", "), "\n")
cat("Class_small categories (", length(class_small_levels), "):", paste(class_small_levels, collapse=", "), "\n")

# ============================================================================
# 2. 数据准备
# ============================================================================

# Extract principal component data
PC1 <- data$PC1
PC2 <- data$PC2  
#PC3 <- data$PC3

# Create a complete data frame
frame <- data.frame(
  ID = data$ID,
  PC1 = PC1,
  PC2 = PC2,
  #PC3 = PC3, #! 如果没有PC3数据，可以注释掉这一行
  Class_big = data$Class_big,
  Class_small = data$Class_small
)

cat("\nData frame created with", nrow(frame), "samples\n")

# ============================================================================
# 3. Color and shape mapping
# ============================================================================


# ========== New: read color mapping from color.csv ===========
color_map_file <- color_file
color_df <- read.csv(color_map_file, header = TRUE, stringsAsFactors = FALSE)
color_map_from_file <- setNames(color_df$color, color_df$Class_big)
# Automatic palette for colors
auto_palette <- brewer.pal(max(8, length(Class_big_levels)), "Set2")
if (length(Class_big_levels) > length(auto_palette)) {
    auto_palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(Class_big_levels))
}
auto_palette_idx <- 1

color_mapping <- c()
for (big_class in Class_big_levels) {
  if (!is.na(color_map_from_file[big_class])) {
    color_mapping[big_class] <- color_map_from_file[big_class]
  } else {
    # Automatically assign color
    color_mapping[big_class] <- auto_palette[auto_palette_idx]
    auto_palette_idx <- auto_palette_idx + 1
    cat(sprintf("No predefined color for category %s, assigning automatically\n", big_class))
  }
}

# Define point shapes (supported by base R plotting)
shape_palette <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                  16, 17, 18, 19, 20, 21, 22, 23, 24, 25)

# Assign shapes to Class_small within each Class_big
create_shape_mapping <- function(data, big_levels, shape_values) {
  shape_map <- c()
  
  for (big_class in big_levels) {
    small_classes <- unique(data$Class_small[data$Class_big == big_class])
    small_classes <- sort(small_classes)  # Sort to ensure consistency
    
    # Restart from the first shape for each big category
    shape_idx <- 1
    
    for (small_class in small_classes) {
      shape_map[small_class] <- shape_values[shape_idx]
      shape_idx <- shape_idx + 1
      if (shape_idx > length(shape_values)) shape_idx <- 1
    }
  }
  return(shape_map)
}

shape_mapping <- create_shape_mapping(data, Class_big_levels, shape_palette)

# Print mapping information
cat("\n=== Color mapping ===\n")
for (i in seq_along(color_mapping)) {
  cat(sprintf("%-12s: %s\n", names(color_mapping)[i], color_mapping[i]))
}

cat("\n=== Shape mapping ===\n")
for (i in seq_along(shape_mapping)) {
  cat(sprintf("%-20s: %2d\n", names(shape_mapping)[i], shape_mapping[i]))
}

# ============================================================================
# 4. Plot parameter settings
# ============================================================================


# Define plotting parameters
plot_params <- list(
  # 图形尺寸
  pdf_width = 35, #* PDF width
  pdf_height = 35, #* PDF height
  
  # Axis parameters
  axis_text_size = 2,      # Axis tick text size
  axis_label_size = 3,     # Axis label size
  title_size = 3,          # Title size
  
  # Scatter point parameters
  point_size = 3,          # Point size
  point_alpha = 0.8,         # Point transparency (via color)
  point_lwd = 3,             # Point border line width
  
  # Legend parameters
  legend_text_size = 2,    # Legend text size
  legend_point_size = 2,   # Point size in legend
  legend_ncol = 5            # Number of legend columns
)

cat("\n=== Plot parameters ===\n")
cat("Point size:", plot_params$point_size, "\n")
cat("Axis text size:", plot_params$axis_text_size, "\n")
cat("Axis label size:", plot_params$axis_label_size, "\n")

# ============================================================================
# 5. Main plot drawing
# ============================================================================

cat("\nStart drawing PCA scatter plot...\n")

# Create PDF file
#! Note: modify the path according to your actual data location
pdf(OUT_pdf_file, width = plot_params$pdf_width, height = plot_params$pdf_height)

# Set plotting layout: main plot on top, legend below
layout(matrix(c(1, 2)), widths = c(1, 1), heights = c(4, 1))

# === Draw main scatter plot ===
par(mar = c(5, 5, 4, 2))  # Increase margins to accommodate larger text

# Create an empty plotting frame
plot(PC1, PC2, 
     type = "n",                           # Do not draw points, just create the frame
     xlab = "PC1",                         # X-axis label
     ylab = "PC2",                         # Y-axis label
     main = "PCA",                         # Title
     cex.axis = plot_params$axis_text_size,    # Axis tick text size
     cex.lab = plot_params$axis_label_size,    # Axis label size
     cex.main = plot_params$title_size,        # Title size
     mgp = c(3, 1, 0))                     # Axis position adjustment



# Draw scatter points by category

# Draw scatter points by category with configurable border width
draw_points_by_class <- function(frame, big_levels, color_map, shape_map, point_size, point_lwd = 1) {
  points_drawn <- 0
  
  for (big_class in big_levels) {
    small_classes <- unique(frame$Class_small[frame$Class_big == big_class])
    small_classes <- sort(small_classes)
    
    cat("Drawing category", big_class, "with", length(small_classes), "subcategories\n")
    
    for (small_class in small_classes) {
      # Subset data for current category
      current_data <- subset(frame, Class_big == big_class & Class_small == small_class)
      
      if (nrow(current_data) > 0) {
        # Draw scatter points
        points(current_data$PC1, current_data$PC2,
               pch = shape_map[small_class],      # Point shape
               bg = color_map[big_class],         # Fill color
               col = color_map[big_class],        # Border color
               cex = point_size,                  # Point size
               lwd = point_lwd)                   # Border line width
        
        points_drawn <- points_drawn + nrow(current_data)
        cat("  -", small_class, ":", nrow(current_data), "points\n")
      }
    }
  }
  
  cat("Total number of points drawn:", points_drawn, "\n")
}

# Execute drawing of scatter points with border width parameter
draw_points_by_class(frame, Class_big_levels, color_mapping, shape_mapping, plot_params$point_size, plot_params$point_lwd)

# ============================================================================
# 6. Legend drawing
# ============================================================================

# === Draw legend ===
par(mar = c(1, 1, 1, 1))  # Reduce margins for the legend area
plot.new()

# Create legend data
create_legend_data <- function(data, big_levels, color_map, shape_map) {
  legend_items <- list(
    labels = c(),
    colors = c(),
    shapes = c(),
    fonts = c()
  )
  
  for (big_class in big_levels) {
    # Add big category title
    legend_items$labels <- c(legend_items$labels, big_class)
    legend_items$colors <- c(legend_items$colors, NA)
    legend_items$shapes <- c(legend_items$shapes, NA)
    legend_items$fonts <- c(legend_items$fonts, 2)  # Bold font
    
    # Add subcategories
    small_classes <- unique(data$Class_small[data$Class_big == big_class])
    small_classes <- sort(small_classes)
    
    for (small_class in small_classes) {
      legend_items$labels <- c(legend_items$labels, paste("  ", small_class))
      legend_items$colors <- c(legend_items$colors, color_map[big_class])
      legend_items$shapes <- c(legend_items$shapes, shape_map[small_class])
      legend_items$fonts <- c(legend_items$fonts, 1)  # Regular font
    }
  }
  
  return(legend_items)
}

legend_data <- create_legend_data(data, Class_big_levels, color_mapping, shape_mapping)

# Draw legend
legend("center", 
       legend = legend_data$labels,
       pch = legend_data$shapes,
       pt.bg = legend_data$colors,
       col = legend_data$colors,
       ncol = plot_params$legend_ncol,
       cex = plot_params$legend_text_size,
       pt.cex = plot_params$legend_point_size,
       bty = "n",                          # No border
       text.font = legend_data$fonts,
       xpd = TRUE)

# Close PDF device
dev.off()

cat("\nPDF file has been saved as: PCA.pdf\n")

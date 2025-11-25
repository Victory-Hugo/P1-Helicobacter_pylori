#!/usr/bin/Rscript

#======================================================================
# TreeMix result analysis and visualization script
# Instructions:
# This script uses predefined paths and can be executed line-by-line in R
# Works with TreeMix output files that follow the Treemix{m}.{r} naming pattern
#======================================================================

# Helper to determine the script directory for relative paths
get_script_dir <- function() {
  if (!interactive()) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- "--file="
    script_path <- sub(file_arg, "", args[grep(file_arg, args)])
    if (length(script_path) > 0) {
      return(dirname(normalizePath(script_path)))
    }
  }
  return(normalizePath(getwd()))
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

# ======= Path and parameter settings (modify as needed) =======

# TreeMix output directory
treemix_dir <- file.path(project_root, "output", "China_only")

# Population list file path
pop.uniq <- file.path(project_root, "conf", "pop_China_only.csv")

# TreeMix output file prefix
file_prefix <- "Treemix"

# Results output directory
results_dir <- file.path(treemix_dir, "analysis")

# Create results directory if missing
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Switch to TreeMix output directory
setwd(treemix_dir)
cat(sprintf("Working directory set to: %s\n", getwd()))

# Load TreeMix plotting functions
source(file.path(project_root, "func", "plotting_funcs.R"))

# Verify that the population list file exists
if (!file.exists(pop.uniq)) {
  cat("Error: Population file not found:", pop.uniq, "\n")
  cat("Please provide a valid population name list file\n")
  stop("Execution stopped because of the missing population file")
}

# ======= File analysis =======
cat("Inspecting output files...\n")
all_files <- list.files(path = treemix_dir, pattern = "covse.gz")

if (length(all_files) == 0) {
  stop(sprintf("No TreeMix output files (*.covse.gz) found in %s. Check the directory or filenames.", treemix_dir))
}

# Extract the different migration edge counts
# Adjust regex to properly match multi-digit migration edge counts
m_values <- unique(gsub(paste0("^", file_prefix, "([0-9]+)\\.[0-9]+\\.covse\\.gz$"), "\\1", all_files))
m_values <- sort(as.numeric(m_values))
cat(sprintf("Detected migration edge counts: %s\n", paste(m_values, collapse = ", ")))

# Determine replicate count for each migration edge number
rep_counts <- sapply(m_values, function(m) {
  files <- list.files(path = treemix_dir, pattern = paste0("^", file_prefix, m, "\\.[0-9]+\\.covse\\.gz$"))
  rep_nums <- as.numeric(gsub(paste0("^", file_prefix, m, "\\.([0-9]+)\\.covse\\.gz$"), "\\1", files))
  max_rep <- max(rep_nums)
  return(max_rep)
})
names(rep_counts) <- m_values

cat("Migration edge counts and replicate numbers:\n")
for (m in m_values) {
  cat(sprintf("  m=%s: %d replicates\n", m, rep_counts[m]))
}

max_m <- max(as.numeric(m_values))
cat(sprintf("Maximum migration edge count: m=%d\n", max_m))

#======================================================================
# Part 1: Generate tree plots for each migration edge count
#======================================================================

# This section temporarily uses the first replicate; the best replicate is selected later
# The best replicate will be determined when computing variance explained
cat("Plotting migration edge trees (using the first replicate)...\n")
for (m in m_values) {
  rep <- 1  # Default to the first replicate; later replaced with the best replicate

  prefix <- sprintf("%s%s.%d", file_prefix, m, rep)
  pdf_file <- file.path(results_dir, paste0("treemix.m", m, ".pdf"))
  cat(sprintf("  Creating figure: %s (using %s)\n", basename(pdf_file), prefix))

  pdf(file = pdf_file, width = 12, height = 8)
  plot_tree(prefix)
  dev.off()
}

#======================================================================
# Part 2: Select the best replicate and compute variance explained
#======================================================================

# Extract likelihood values from .llik files
extract_likelihood <- function(llik_file) {
  if (!file.exists(llik_file)) {
    return(NA)
  }

  lines <- readLines(llik_file)
  # Identify the "Exiting ln(likelihood)" line, which typically contains the final likelihood
  exit_line <- grep("Exiting ln\\(likelihood\\)", lines, value = TRUE)

  if (length(exit_line) > 0) {
    # Extract the likelihood value (the last number)
    likelihood <- as.numeric(gsub(".*: ([0-9.-]+).*", "\\1", exit_line))
    return(likelihood)
  }

  return(NA)
}

# Identify the best replicate per migration edge count
cat("Identifying the replicate with the highest likelihood for each migration edge count...\n")
best_reps <- list()

# Ensure variables are pre-initialized
for (m in m_values) {
  # Safely initialize the best replicate
  best_rep <- 1  # Default to the first replicate

  # Attempt to collect likelihoods for all replicates
  tryCatch({
    rep_count <- as.numeric(rep_counts[as.character(m)])
    if (is.na(rep_count) || rep_count < 1) {
      cat(sprintf("  Warning: m=%s has no valid replicate count; defaulting to rep=1\n", m))
    } else {
      likelihoods <- numeric(rep_count)

      # Collect all available likelihood values
      for (r in 1:rep_count) {
        llik_file <- file.path(treemix_dir, sprintf("%s%s.%d.llik", file_prefix, m, r))
        if (file.exists(llik_file)) {
          likelihoods[r] <- extract_likelihood(llik_file)
        } else {
          cat(sprintf("  Warning: m=%s rep=%d lacks an llik file\n", m, r))
          likelihoods[r] <- NA
        }
      }

      # Identify the replicate with the highest likelihood
      if (all(is.na(likelihoods))) {
        cat(sprintf("  Warning: m=%s has no available likelihood values; defaulting to rep=1\n", m))
      } else {
        max_idx <- which.max(likelihoods)
        if (length(max_idx) > 0 && !is.na(likelihoods[max_idx])) {
          best_rep <- max_idx
          cat(sprintf("  m=%s: best replicate is rep=%d (likelihood: %.2f)\n",
                      m, best_rep, likelihoods[best_rep]))
        } else {
          cat(sprintf("  Warning: m=%s likelihoods are all NA; defaulting to rep=1\n", m))
        }
      }
    }
  }, error = function(e) {
    cat(sprintf("  Error: Failed to process m=%s: %s\n", m, e$message))
    best_rep <<- 1  # Use the first replicate if an error occurs
  })

  best_reps[[as.character(m)]] <- best_rep
  cat(sprintf("  Selected replicate rep=%d for m=%s\n", best_rep, m))
}

cat("\nComputing variance explained for each migration edge count...\n")

# Compute the variance explained using the best replicate per migration edge count
m_results <- numeric(length(m_values))
names(m_results) <- m_values
best_prefixes <- character(length(m_values))
names(best_prefixes) <- m_values

for (i in 1:length(m_values)) {
  m <- m_values[i]
  rep <- best_reps[[as.character(m)]]

  prefix <- sprintf("%s%s.%d", file_prefix, m, rep)
  best_prefixes[i] <- prefix
  m_results[i] <- get_f(prefix)
  cat(sprintf("  m=%s (rep=%d) variance explained: %.4f%%\n", m, rep, m_results[i] * 100))
}

# Save variance explained results
results_table <- data.frame(
  migration_edges = m_values,
  best_rep = sapply(m_values, function(m) best_reps[[as.character(m)]]),
  best_prefix = best_prefixes,
  variance_explained = m_results * 100
)
csv_file <- file.path(results_dir, "variance_explained.csv")
write.table(results_table, csv_file,
            quote = FALSE, sep = ",", row.names = FALSE)
cat(sprintf("Variance explained table saved to %s\n", csv_file))

# Determine the smallest migration edge count that reaches 99.8%
threshold <- 99.8
threshold_m <- min(m_values[m_results * 100 >= threshold], Inf)
if (is.infinite(threshold_m)) {
  cat(sprintf("Warning: No migration edge count reaches %.1f%% variance explained\n", threshold))
  cat(sprintf("Maximum variance explained: %.4f%% (m=%s)\n", max(m_results) * 100, m_values[which.max(m_results)]))
} else {
  cat(sprintf("Smallest migration edge count reaching %.1f%% variance explained: m=%d\n", threshold, threshold_m))

  # Copy the best model file
  best_rep <- best_reps[[as.character(threshold_m)]]
  best_tree_file <- sprintf("%s%s.%d.treeout.gz", file_prefix, threshold_m, best_rep)
  if (file.exists(file.path(treemix_dir, best_tree_file))) {
    file.copy(
      from = file.path(treemix_dir, best_tree_file),
      to = file.path(results_dir, "best_model.treeout.gz"),
      overwrite = TRUE
    )
    cat(sprintf("Copied best model file to: %s\n", file.path(results_dir, "best_model.treeout.gz")))

    # Regenerate the tree plot for the best model
    best_prefix <- sprintf("%s%s.%d", file_prefix, threshold_m, best_rep)
    pdf_file <- file.path(results_dir, paste0("best_treemix_m", threshold_m, "_rep", best_rep, ".pdf"))
    cat(sprintf("Generating best model figure: %s\n", basename(pdf_file)))
    pdf(file = pdf_file, width = 12, height = 8)
    plot_tree(best_prefix)
    dev.off()
  }
}

# Plot migration edge count versus variance explained
cat("Plotting migration edge selection chart...\n")
pdf_file <- file.path(results_dir, "migration.edge.choice.pdf")
pdf(file = pdf_file, width = 10, height = 7)
plot(as.numeric(m_values), m_results * 100,
     pch = 19, cex = 1.2, col = "blue", type = "b",
     xlab = "Migration edge count", ylab = "Variance explained (%)",
     xaxt = "n",
     main = "TreeMix variance explained")

# Custom x-axis labels
axis(1, at = as.numeric(m_values), labels = m_values)

# 99.8% reference line
abline(h = threshold, col = "red", lty = 2)
text(as.numeric(m_values)[2], threshold + 0.5,
     sprintf("%.1f%% threshold", threshold), col = "red")

# Highlight the smallest m that reaches the threshold
if (!is.infinite(threshold_m)) {
  points(threshold_m, m_results[as.character(threshold_m)] * 100,
         pch = 19, col = "red", cex = 2)
  text(threshold_m, m_results[as.character(threshold_m)] * 100 - 1.5,
       sprintf("m=%d", threshold_m), col = "red")
}

# Add variance explained labels
for (i in 1:length(m_values)) {
  text(as.numeric(m_values)[i], m_results[i] * 100 + 0.7,
       sprintf("%.2f%%", m_results[i] * 100), cex = 0.8)
}

dev.off()
cat(sprintf("Figure saved to: %s\n", pdf_file))

#======================================================================
# Part 3: Residual analysis
#======================================================================

cat("Generating residual heatmaps...\n")
# Verify that the population file matches TreeMix output files
check_pop_names <- function(prefix, pop_file) {
  cov_file <- paste0(prefix, ".cov.gz")
  if (!file.exists(cov_file)) {
    cat(sprintf("Warning: File %s does not exist\n", cov_file))
    return(FALSE)
  }

  tryCatch({
    pop_names <- readLines(pop_file)
    cat(sprintf("Population file contains %d population names\n", length(pop_names)))
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("Error while reading population file: %s\n", e$message))
    return(FALSE)
  })
}

# Draw residual plots for each migration edge count (best replicate + error handling)
cat("Generating residual heatmaps...\n")
for (m in m_values) {
  tryCatch({
    rep <- best_reps[[as.character(m)]]

    prefix <- sprintf("%s%s.%d", file_prefix, m, rep)
    pdf_file <- file.path(results_dir, sprintf("treemix.residuals.m%s.pdf", m))
    cat(sprintf("  Attempting to create figure: %s (using %s)\n", basename(pdf_file), prefix))

    if (check_pop_names(prefix, pop.uniq)) {
      pdf(file = pdf_file)
      tryCatch({
        plot_resid(prefix, pop.uniq)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, sprintf("Error while plotting residuals:\n%s", e$message), cex = 1.2)
        cat(sprintf("Error: %s\n", e$message))
      })
      dev.off()

      if (!is.infinite(threshold_m) && m == threshold_m) {
        hi_res_file <- file.path(results_dir, "best_model_residuals.pdf")
        cat(sprintf("  Generating best model residual plot: %s\n", basename(hi_res_file)))
        pdf(file = hi_res_file, width = 12, height = 10)
        tryCatch({
          plot_resid(prefix, pop.uniq)
        }, error = function(e) {
          plot.new()
          text(0.5, 0.5, sprintf("Error while plotting best-model residuals:\n%s", e$message), cex = 1.2)
        })
        dev.off()
      }
    }
  }, error = function(e) {
    cat(sprintf("Error while processing m=%s: %s\n", m, e$message))
  })
}

# Create summary report
report_file <- file.path(results_dir, "treemix_analysis_summary.txt")
cat("Creating analysis summary report...\n")

report_con <- file(report_file, "w")
cat("=========================================\n", file = report_con)
cat("TreeMix variance explained analysis summary\n", file = report_con)
cat("=========================================\n\n", file = report_con)
cat(sprintf("Analysis date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = report_con)
cat("Analyzed migration edge counts:\n", file = report_con)
for (i in 1:length(m_values)) {
  cat(sprintf("  m=%s: %.4f%%\n", m_values[i], m_results[i] * 100), file = report_con)
}
cat("\n", file = report_con)

if (!is.infinite(threshold_m)) {
  best_rep <- best_reps[[as.character(threshold_m)]]
  cat(sprintf("Smallest migration edge count reaching %.1f%% variance explained: m=%d\n", threshold, threshold_m), file = report_con)
  cat(sprintf("Best model: %s%s.%d\n", file_prefix, threshold_m, best_rep), file = report_con)
  cat(sprintf("Variance explained: %.4f%%\n", m_results[as.character(threshold_m)] * 100), file = report_con)
} else {
  best_m <- m_values[which.max(m_results)]
  best_rep <- best_reps[[as.character(best_m)]]
  cat(sprintf("Warning: No migration edge count reaches %.1f%% variance explained\n", threshold), file = report_con)
  cat(sprintf("Maximum variance explained: %.4f%% (m=%s, rep=%d)\n",
              max(m_results) * 100, best_m, best_rep), file = report_con)
}
close(report_con)

cat("All analyses completed!\n")
cat(sprintf("Results summarize %d migration edge counts\n", length(m_values)))
cat(sprintf("All results saved in: %s\n", results_dir))

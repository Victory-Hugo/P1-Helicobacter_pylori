#!/usr/bin/env Rscript

## 1) Load packages
suppressMessages({
  if (!requireNamespace("xml2", quietly=TRUE)) install.packages("xml2")
  if (!requireNamespace("ape",   quietly=TRUE)) install.packages("ape")
  if (!requireNamespace("treeio",quietly=TRUE)) install.packages("treeio")
})
library(xml2)
library(ape)
library(treeio)

## 2) Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: convert_tree.R <input tree.xml path> <output .nwk path>\n")
  q(status = 1)
}
treefile <- args[1]
outfile  <- args[2]

## 3) Read XML and extract <Tree> node text
xml_doc   <- read_xml(treefile)
tree_text <- xml_text(xml_find_all(xml_doc, ".//Tree"))

if (length(tree_text) == 0) {
  stop("No <Tree> node found in the XML file. Check the input file.")
}

## 4) Parse to phylo object and save Newick file
phy <- read.tree(text = tree_text)
write.tree(phy, file = outfile)

cat("Converted the tree in", treefile, "to Newick format at", outfile, "\n")

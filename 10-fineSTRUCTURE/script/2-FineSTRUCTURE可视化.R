###########################################
#### FineSTRUCTURE analysis script (optimized example) ####
###########################################

## Notes:
#todo 1. Install FineSTRUCTURE v2 and ensure the fs executable is on PATH or referenced explicitly.
#todo 2. Install required R packages such as "gplots" and "XML".
#todo 3. Adjust paths and filenames based on your dataset.

#################################
#### 1. Load R packages and helper functions ####
#################################
rm(list=ls())
# conda install -c conda-forge r-gplots
# install.packages("ape")
# install.packages("XML")
library("gplots") 

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
output_dir <- file.path(project_root, "output")
python_dir <- file.path(project_root, "python")
heatmap_script <- file.path(python_dir, "heatmap.3.R")
matrix_lib <- file.path(python_dir, "Matrix_FinestructureLibrary.R")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

source(matrix_lib)
## The referenced FinestructureLibrary.R loads the helper functions and dependencies needed for this analysis


##########################################
#### 2. Define FineSTRUCTURE output paths ####
##########################################
## The following three files come from the main FineSTRUCTURE workflow
chunkfile <- "WGSHP_linked.chunkcounts.out"
mcmcfile  <- "WGSHP_linked_mcmc.xml"
treefile  <- "WGSHP_linked_tree.xml"

########################################################
#### 3. Generate additional FineSTRUCTURE outputs ####
########################################################
## (1) population-by-population chunkcount file
mappopchunkfile <- "WGSHP_linked.mutationprobs.out"
## Run fs with -X -Y -e X2

command_X2 <- paste("fs fs -X -Y -e X2", chunkfile, treefile, mappopchunkfile)
system(command_X2)

## (2) Pairwise coincidence file (mean co-clustering across the MCMC)
meancoincidencefile <- "WGSHP_linked.meancoincidence.csv"
## Run fs with -X -Y -e meancoincidence
command_meancoin <- paste("fs fs -X -Y -e meancoincidence", chunkfile, mcmcfile, meancoincidencefile)
system(command_meancoin)

##########################################
#### 4. Read chunkcounts, MCMC, and tree data ####
##########################################
## 4.1 Read chunkcounts (coancestry matrix)
dataraw <- as.matrix(read.table(chunkfile, 
                                row.names = 1, 
                                header = TRUE, 
                                skip = 1))

## 4.2 Parse the MCMC XML file
library("XML")
mcmcxml  <- xmlTreeParse(mcmcfile)
mcmcdata <- as.data.frame.myres(mcmcxml)

## 4.3 Parse and extract the tree XML
treexml <- xmlTreeParse(treefile)
ttree   <- extractTree(treexml)     # Convert to ape::phylo
ttree$node.label <- NULL            # Remove internal labels to avoid clutter
tdend   <- myapetodend(ttree, 
                       factor = 1)  # Convert phylo object to a dendrogram

###########################################
#### 5. Extract and process FineSTRUCTURE clusters ####
###########################################
## Obtain MAP (Maximum A Posteriori) group assignments
mapstate     <- extractValue(treexml, "Pop")  # Cluster assignments
mapstatelist <- popAsList(mapstate)           # List individuals per cluster

## Create cluster names
popnames      <- lapply(mapstatelist, NameSummary)      # Reversible format (lossless)
popnamesplot  <- lapply(mapstatelist, NameMoreSummary)  # Readable visualization labels
names(popnames)     <- popnamesplot
names(popnamesplot) <- popnamesplot

## Build dendrograms and adjust midpoint membership
popdend      <- makemydend(tdend, mapstatelist)
popdend      <- fixMidpointMembers(popdend)
popdendclear <- makemydend(tdend, mapstatelist, "NameMoreSummary")
popdendclear <- fixMidpointMembers(popdendclear)

##################################################
#### 6. Read and process the pairwise coincidence matrix ####
##################################################
## Reorder the coincidence matrix to match tree order
fullorder       <- labels(tdend) 
mcmcmatrixraw   <- as.matrix(read.csv(meancoincidencefile, 
                                      row.names = 1))
mcmcmatrix      <- mcmcmatrixraw[fullorder, fullorder]

## Retrieve the MAP-state binary matrix (groupingAsMatrix)
mapstatematrix  <- groupingAsMatrix(mapstatelist)[fullorder, fullorder]

######################################################
#### 7. Process and store the coancestry matrix ####
######################################################
datamatrix <- dataraw[fullorder, fullorder]

## Export the coancestry matrix to CSV for downstream ID replacement
datamatrix_path <- file.path(output_dir, "datamatrix.csv")
write.csv(datamatrix, datamatrix_path)

## Note: external tools (e.g., Excel VLOOKUP) can be used to replace INDXXX with actual sample IDs

############################################################
#!### 8. Reload the matrix after ID replacement ####
#! Example CLI: awk -v FS=',' '{gsub("\"", "", $1); print $1}' "<datamatrix_path>" > sample_list.csv
############################################################
InputFileData  <- datamatrix_path
datamatrix_R   <- read.table(file = InputFileData, 
                             header = TRUE, 
                             sep = ",", 
                             check.names = FALSE)

## Restore first column as row names
names_pop            <- datamatrix_R[, 1]
datamatrix_R[, 1]    <- NULL
rownames(datamatrix_R) <- names_pop
datamatrix_R         <- as.matrix(datamatrix_R)

##################################################
#### 9. Define custom color functions ####
##################################################
MakeCustomColors <- function(numColors = 300, useDarkEnd = FALSE) {
  ## Define a base palette
  baseColors <- c("#003936", "#31646C", "#4E9280", "#96B89B", 
                  "#DCDFD2", "#D49C87", "#B86265", "#50184E")
  
  ## Optionally force a dark ending color
  if (useDarkEnd) {
    baseColors[length(baseColors)] <- "#333333"  # Replace last color with dark gray
  }

  ## Build gradient function
  colorFunction <- colorRampPalette(baseColors)
  
  ## Generate the requested number of colors
  return(colorFunction(numColors))
}

## Example palettes
some.colors    <- MakeCustomColors()
some.colorsEnd <- MakeCustomColors(numColors = 100, TRUE)  # Palette with dark gray ending

## Cap extreme values at tmatmax for better visualization
tmatmax       <- 200
tmpmat        <- datamatrix_R
tmpmat[tmpmat > tmatmax] <- tmatmax



######################################################
#### 10. Load sample annotations and define colors ####
######################################################
## Assign colors for rows / columns based on population/group metadata
cols      <- rep('#828282', ncol(datamatrix_R))  # Default color
cols_rows <- rep('#828282', nrow(datamatrix_R))
#! Update the following as needed
list_color <- read.table(file = "list.csv", 
                         header = TRUE, 
                         sep = ",")

## Assign colors by POP
Population <- rep('#CEF9FF', nrow(list_color))

Population[list_color$POP %in% "Han"]                     <- '#FDD262'
Population[list_color$POP %in% "hpAfrica1"]               <- '#E2D200'
Population[list_color$POP %in% "hpAfrica2"]               <- '#CB2314'
Population[list_color$POP %in% "hpAsia2"]                 <- '#79402E'
Population[list_color$POP %in% "hspEAsia"]                <- '#0B775E'
Population[list_color$POP %in% "hpEurope"]                <- '#DC863B'
Population[list_color$POP %in% "hpEuropeSahul"]           <- '#FDDDA0'
Population[list_color$POP %in% "hpNEAfrica"]              <- '#C27D38'
Population[list_color$POP %in% "hpNorthAsia"]             <- '#9B110E'
Population[list_color$POP %in% "hpSahul"]                 <- '#C52E19'
Population[list_color$POP %in% "hspAfrica1MiscAmerica"]   <- '#E79805'
Population[list_color$POP %in% "hspAfrica1NAmerica"]      <- '#D3DDDC'
Population[list_color$POP %in% "hspAfrica1SAfrica"]       <- '#E6A0C4'
Population[list_color$POP %in% "hspAfrica1WAfrica"]       <- '#ABDDDE'
Population[list_color$POP %in% "hspEurasia"]               <- '#ECCBAE'
Population[list_color$POP %in% "hspIndigenousAmerica"]    <- '#E1AF00'
Population[list_color$POP %in% "hspNEurope"]               <- '#F8AFA8'
Population[list_color$POP %in% "hspSWEurope"]              <- '#CCBA72'
Population[list_color$POP %in% "hspUral"]                  <- '#0E2C68'
Population[list_color$POP %in% "T"]                        <- '#000000'

## Assign colors by DAPC clusters
colours_list_DAPC <- rep('#CEF9FF', nrow(list_color))
colours_list_DAPC[list_color$DAPC %in% "East_Asia"] <- '#008041'
colours_list_DAPC[list_color$DAPC %in% "Southeast_Asia"] <- '#0081cf'
colours_list_DAPC[list_color$DAPC %in% "South_Asia"] <- '#1951af'
colours_list_DAPC[list_color$DAPC %in% "Other"] <- '#000000'
colours_list_DAPC[list_color$DAPC %in% "Oceania"] <- '#aB3D4F'
colours_list_DAPC[list_color$DAPC %in% "Europe"] <- '#d11fae'
colours_list_DAPC[list_color$DAPC %in% "Central_Asia"] <- '#C1D2D2'
colours_list_DAPC[list_color$DAPC %in% "Papua_New_Guinea"] <- '#CB3D4F'
colours_list_DAPC[list_color$DAPC %in% "America"] <- '#947145'
colours_list_DAPC[list_color$DAPC %in% "South_Africa"] <- '#7F7F7F'
colours_list_DAPC[list_color$DAPC %in% "North_America"] <- '#DD8D29'
colours_list_DAPC[list_color$DAPC %in% "South_America"] <- '#FAD510'
colours_list_DAPC[list_color$DAPC %in% "West_Africa"] <- '#FAEFD1'
colours_list_DAPC[list_color$DAPC %in% "Central_Africa"]  <- '#9AA83A'
colours_list_DAPC[list_color$DAPC %in% "West_Asia"]  <- '#00668c'
colours_list_DAPC[list_color$DAPC %in% "Europe-Asia"]  <- '#81A88D'
colours_list_DAPC[list_color$DAPC %in% "North_Africa"]  <- '#000000'
# # Assign colors by additional grouping labels
# colours_list_Classification <- rep('#7F7F7F', nrow(list_color))
# colours_list_Classification[list_color$Classification %in% "Other"] <- '#7F7F7F'
# colours_list_Classification[list_color$Classification %in% "Russia"] <- '#81A88D'
# colours_list_Classification[list_color$Classification %in% "Korea"] <- '#FAD510'
# colours_list_Classification[list_color$Classification %in% "Mongolia"] <- '#DD8D29'
# colours_list_Classification[list_color$Classification %in% "China"] <- '#de283b'
# colours_list_Classification[list_color$Classification %in% "Australia"] <- '#00668c'
# colours_list_Classification[list_color$Classification %in% "Indonesia"] <- '#1951af'
# colours_list_Classification[list_color$Classification %in% "India"] <- '#71c4ef'
# colours_list_Classification[list_color$Classification %in% "Thailand"] <- '#afffff'
# colours_list_Classification[list_color$Classification %in% "Myanmar"] <- '#25b1bf'
# colours_list_Classification[list_color$Classification %in% "Malaysia"] <- '#008041'
# colours_list_Classification[list_color$Classification %in% "Japan"] <- '#ffaea0'
# colours_list_Classification[list_color$Classification %in% "Nepal"] <- '#228B22'
# colours_list_Classification[list_color$Classification %in% "Vietnam"] <- '#FAEFD1'
# colours_list_Classification[list_color$Classification %in% "Bhutan"] <- '#3A9AB2'

# Assign colors by coarse classification
colours_list_Classification <- rep('#7F7F7F', nrow(list_color))
colours_list_Classification[list_color$Classification %in% "South"]  <- '#1951af'
colours_list_Classification[list_color$Classification %in% "North"] <- '#de283b'


colours_list_Classification2 <- rep('#7F7F7F', nrow(list_color))
colours_list_Classification2[list_color$Classification2 %in% "Other"] <- '#000000'
colours_list_Classification2[list_color$Classification2 %in% "Inner_Mongolia"] <- '#CB3D4F'
colours_list_Classification2[list_color$Classification2 %in% "Heilongjiang"] <- '#d11fae'
colours_list_Classification2[list_color$Classification2 %in% "Ningxia"] <- '#FF3D3D'
colours_list_Classification2[list_color$Classification2 %in% "Beijing"] <- '#de283b'
colours_list_Classification2[list_color$Classification2 %in% "Sichuan"] <- '#00668c'
colours_list_Classification2[list_color$Classification2 %in% "Xizang"] <- '#1951af'
colours_list_Classification2[list_color$Classification2 %in% "Guangxi"] <- '#71c4ef'
colours_list_Classification2[list_color$Classification2 %in% "Yunnan"] <- '#afffff'
colours_list_Classification2[list_color$Classification2 %in% "Guizhou"] <- '#25b1bf'
colours_list_Classification2[list_color$Classification2 %in% "Fujian"] <- '#5B2C6F'
colours_list_Classification2[list_color$Classification2 %in% "Zhejiang"] <- '#3F51B5'
colours_list_Classification2[list_color$Classification2 %in% "Taiwan"] <- '#757de8'
colours_list_Classification2[list_color$Classification2 %in% "Shenzhen"] <- '#0E549B'
colours_list_Classification2[list_color$Classification2 %in% "Hongkong"] <- '#9AA83A'
colours_list_Classification2[list_color$Classification2 %in% "Shandong"] <- '#6E211F'
colours_list_Classification2[list_color$Classification2 %in% "Shanghai"] <- '#52A400'
colours_list_Classification2[list_color$Classification2 %in% "Hunan"] <- '#0E73CF'
colours_list_Classification2[list_color$Classification2 %in% "Chengdu"] <- '#00668c'
colours_list_Classification2[list_color$Classification2 %in% "Gansu"] <- '#DC5C00'
colours_list_Classification2[list_color$Classification2 %in% "Shaanxi"] <- '#FF4E38'
colours_list_Classification2[list_color$Classification2 %in% "Taipei"] <- '#757de8'


# Build row/column annotation matrices for heatmap side colors
#! Customize as required; uncomment additional labels when needed
rlab=as.matrix(t(cbind(colours_list_ClassificatioWGSn)))
clab=as.matrix((cbind(colours_list_Classification)))
rlab=as.matrix(t(cbind(colours_list_Classification,colours_list_Classification2)))
clab=as.matrix((cbind(colours_list_Classification,colours_list_Classification2)))
rlab=as.matrix(t(cbind(Population,colours_list_DAPC)))
clab=as.matrix((cbind(Population,colours_list_DAPC)))

# rlab=as.matrix(t(cbind(Population,colours_list_DAPC,colours_list_Classification)))
# clab=as.matrix((cbind(Population,colours_list_DAPC,colours_list_Classification)))

# rlab=as.matrix(t(cbind(Population,colours_list_DAPC,colours_list_Classification,colours_list_Classification2)))
# clab=as.matrix((cbind(Population,colours_list_DAPC,colours_list_Classification,colours_list_Classification2)))
###############################
#### 11. Plot heatmaps and export ####
###############################
# Use the custom heatmap.3 function
dev.off()
pdf_path <- file.path(output_dir, "Fs_coancestry_heatmap.pdf")
pdf(pdf_path,
    height = 100, 
    width  = 100)


source(heatmap_script)

p <- heatmap.3(tmpmat,
               scale            = "none",
               Rowv             = tdend,
               Colv             = tdend,
               dendrogram       = "both",
               trace            = "none",
               key              = TRUE,
               keysize          = 0.5,
               cexRow           = 0.3,
               cexCol           = 0.3,
               col              = some.colorsEnd,  # Palette with dark gray ending
               RowSideColorsSize= 1,
               ColSideColorsSize= 1,
               ColSideColors    = clab,
               RowSideColors    = rlab
)
dev.off()

#! If the PDF is too large, create a bitmap alternative
dev.off()
# Use tiff() instead of pdf() when a raster figure is preferred
tiff_path <- file.path(output_dir, "Fs_coancestry_heatmap.tif")
tiff(tiff_path,
     width = 20,      
     height = 20,     
     units = "in",    
     res = 1000,       
     compression = "lzw")

# Optionally adjust margins
par(mar = c(4,4,2,2))

# Load the custom function (same script reused above)
source(heatmap_script)

# Plot
p <- heatmap.3(tmpmat,
               scale            = "none",
               Rowv             = tdend,
               Colv             = tdend,
               dendrogram       = "both",
               trace            = "none",
               key              = TRUE,
               keysize          = 0.5,
               cexRow           = 0.3,
               cexCol           = 0.3,
               col              = some.colorsEnd,
               RowSideColorsSize= 1,
               ColSideColorsSize= 1,
               ColSideColors    = clab,
               RowSideColors    = rlab
)
dev.off()

#setwd("C:/Users/user/Desktop")
setwd("K:/Workbench/Group_support/Katie")
# Name of the Excel file in the given working directory
infileName <- "180319 Lipid Metabolite Data.xlsx"
sheet <- 2
# Name to give to the heatmap that will be created in the given working directory
heatmapName <- c("180319 Lipid Metabolite Data log10.jpeg")
# Rows you are interested in (comment it out if you want to show them all in the heatmap)
#rowSelection <- c()
# Columns you are interested in (comment it out if you want to show them all in the heatmap)
#colSelection <- c()
# Do row normalization?
normalizeRows <- FALSE
# Do you have at least one sample with two or more replicates?
replicates <- TRUE
# Apply Log10 to the data INSTEAD of normalizing?
applyLog10 <- TRUE
# Do column normalization?
normalizeCols <- FALSE
# Number of categories in the Excel file
numCategories <- 0
# Extract categories from lipid names? E.g. PE(18:0/16:1) has category PE
extractCategories <- FALSE
categoryTitles <- c("Phospholipids")
# Set the following varible to TRUE if you want to display the categories in the heatmap,
# change it to FALSE otherwise
showCategories <- TRUE
# Minimum difference between two values to change the color assigned in the heatmap
bucketSize <- 0.01

if (applyLog10 && normalizeCols) {
    warning("The data will be transformed by log10 function but not normalised by columns.")
}

# Install the required packages not downloaded in user"s R library
packages <- c("tools", "readxl", "tidyverse", "data.table", "pheatmap", "RColorBrewer")
missingPackages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(missingPackages)) {
    install.packages(missingPackages)
}
# Load the required packages, suppressing the output
invisible(lapply(packages, require, character.only = TRUE))

# Implementation of method  is.nan(data) for type 'list'
is.nan.data.frame <- function(x) {
    do.call(cbind, lapply(x, is.nan))
}

# Read in the data file
if (tolower(file_ext(infileName)) %in% c("xls", "xlsx")) {
    data <- read_excel(infileName, sheet=sheet)
} else {
    data <- read.csv(infileName)
}
srcColNames <- colnames(data)
#colnames(data) <- enc2native(colnames(data))
colnames(data) <- gsub("<U\\+[0-9A-Z]+>", "*", enc2native(colnames(data)))

# Detect columns that have the same name as another in the original file
dupCols <- which(grepl("__[0-9]+$", colnames(data)))
if (length(dupCols) > 0) {
    print("The following columns have been renamed to avoid duplications:")
    for (index in dupCols) {
        colName <- colnames(data)[index]
        print(paste0(colName, " (column ", index, ")"))
    }
}

# Set a fixed name for the first column (sample names)
colnames(data)[1] <- c("samples")

# Extract categories from column names
if (extractCategories) {
    categories <- data.frame(t(sapply(strsplit(srcColNames, "\\("), "[[", 1)),
                             stringsAsFactors = FALSE)
    colnames(categories) <- srcColNames
    categories[1, 1] <- categoryTitles
}
# Append the categories included in the data file and remove them from the dataframe
if (numCategories > 0) {
    if (exists("categories")) {
        categories <- rbind(categories, data[seq(1, numCategories),])
    } else {
        categories <- data[seq(1, numCategories),]
    }
    data <- data[-seq(1, numCategories), ]
    data <- mutate_at(data, vars(-samples), as.double)
}
# Update the number of categories if needed
if (extractCategories) {
    numCategories <- numCategories + 1
}

# Keep only the rows the user is interested in
if (exists("rowSelection")) {
    rowSelection <- rowSelection - 1 - numCategories
    data <- data[rowSelection,]
}

# Keep only the columns the user is interested in
if (exists("colSelection")) {
    colSelection <- append(colSelection, 1, after = 0)
    data <- data[, colSelection]
    if (numCategories > 0) {
        categories <- categories[, colSelection]
    }
}

# Replace NA by 0s
data[is.na(data)] <- 0

# Row normalization (based on the sum of each row)
if (normalizeRows) {
    rowOrder <- data$samples
    colOrder <- colnames(data)
    data <- data %>%
        gather(lipids, values, -samples) %>%
        spread(samples, values) %>%
        mutate_if(~ is.numeric(.x), ~ .x / sum(.x, na.rm = TRUE)) %>%
        gather(samples, values, -lipids) %>%
        spread(lipids, values)
    data <- data[match(rowOrder, data$samples),]
    setcolorder(data, colOrder)
    data[is.nan(data)] <- 0
}

# Get the average of the replicates (if there is more than 1)
if (replicates) {
    data <- data %>%
        mutate(samples = trimws(gsub("[0-9]+$", "", samples)))
    rowOrder <- unique(data$samples)
    data <- data %>%
        group_by(samples) %>%
        summarise_all(funs(mean), na.rm = TRUE)
    data <- data[match(rowOrder, data$samples),]
    data[is.nan(data)] <- 0
}

# Column normalization (based on the maximum value of each column)
if (applyLog10) {
    #data[data == 0] <- 1e-50
    data <- data %>%
        mutate_if(~ is.numeric(.x), ~ log10(.x))
} else if (normalizeCols) {
    data <- data %>%
        mutate_if(~ is.numeric(.x), ~ .x / max(.x, na.rm = TRUE))
    data[is.nan(data)] <- 0
}

# Get columns with all 0s (excluding "samples" column)
zeroCols <- colSums(data[-c(1)]) == 0
# Add logic value for "samples" column
zeroCols <- append(zeroCols, FALSE, after = 0)
names(zeroCols)[1] <- "samples"
print("Columns with all zeros:")
print(paste(srcColNames[zeroCols], collapse = ", "))
# Remove columns with all 0s
#data <- data[, !zeroCols]
#if (numCategories > 0) {
#    categories <- categories[, !zeroCols]
#}
zeroCols <- rep(FALSE, ncol(data))

# Generate a matrix from the data (ignoring the first column corresponding to the row names)
dataMatrix <- as.matrix(data[,-1])
# Set the same row and column names for the matrix
rownames(dataMatrix) <- data$samples
colnames(dataMatrix) <- srcColNames[!zeroCols][-c(1)]

if (normalizeCols) {
    breaksList = seq(0, 1, by = bucketSize)
} else {
    minValue <- min(dataMatrix)
    maxValue <- max(dataMatrix)
    if (minValue == -Inf) {
        minValue <- min(dataMatrix[is.finite(dataMatrix)])
        floorValue <- -4
        dataMatrix[!is.finite(dataMatrix)] <- floorValue
    }
    breaksList <- seq(minValue, maxValue, by = bucketSize)
    if (exists("floorValue")) {
        breaksList <- append(breaksList, floorValue, after = 0)
    }
    # If maxValue %% bucketSize != 0, add the 'maxValue' at the end of 'breaksList'
    lastValue <- tail(breaksList, n = 1)
    if (lastValue < maxValue) {
        breaksList <- append(breaksList, lastValue + bucketSize)
    }
}
colorPaletter <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
if ((numCategories > 0) && showCategories) {
    categoryTitles <- as.character(unlist(categories[,1]))
    categories <- categories[,-c(1)]
    # Set the category for each column
    categoryNames <- data.frame(t(categories), row.names=srcColNames[!zeroCols][-c(1)])
    colnames(categoryNames) <- c(categoryTitles)
    # Get a color scheme and give it to the first category
    colorSet <- brewer.pal(12, name="Set3")
    #colorSet <- palette(c("#25D9E2", "#FD928E"))
    catSet <- sort(unique(unlist(categories[1,])))
    # Remove NA category (if there is one)
    catSet <- catSet[!is.na(catSet)]
    names(colorSet) <- catSet
    categoryColors <- list(c(head(colorSet, length(catSet))))

    if (numCategories > 1) {
        schemes <- c("Dark2", "Set1", "Set2", "Accent")
        # Get a color scheme and give one color per category
        for (i in 2:numCategories) {
            colorSet <- brewer.pal(8, name=schemes[i-1])
            catSet <- unique(unlist(categories[i,]))
            # Remove NA category (if there is one)
            catSet <- catSet[!is.na(catSet)]
            names(colorSet) <- catSet
            categoryColors[[i]] <- c(head(colorSet, length(catSet)))
        }
    }
    names(categoryColors) <- c(categoryTitles)

    # Create and save the heatmap figure
    pheatmap(dataMatrix, border_color = NA, cellwidth = 9, cellheight = 9,
             annotation_col = rev(categoryNames), annotation_colors = rev(categoryColors),
             cluster_rows = FALSE, cluster_cols = TRUE, show_colnames = TRUE,
             fontsize = 8, breaks = breaksList, color = colorPaletter,
             filename = heatmapName)#, gaps_row = seq(7, nrow(dataMatrix)-1, by=7))
             #gaps_row = c(3,5,8,11,14))
} else {
    # Create and save the heatmap figure
    pheatmap(dataMatrix, border_color = NA, cellwidth = 9, cellheight = 9,
             cluster_rows = FALSE, cluster_cols = TRUE, show_colnames = TRUE,
             fontsize = 8, breaks = breaksList, color = colorPaletter,
             filename = heatmapName, cutree_col = 5)#, gaps_row = seq(7, nrow(dataMatrix), by=7))
}
# height = 3
# display_numbers = pvMatrix
#   pvalues <- data
#   pvalues <- replace(pvalues, as.matrix(data >= 0.2), "")
#   pvalues <- replace(pvalues, as.matrix(data < 0.2), "*")
#   pvalues <- replace(pvalues, as.matrix(data < 0.1), "**")
#   pvalues <- replace(pvalues, as.matrix(data < 0.05), "***")
#   pvMatrix <- as.matrix(pvalues[, -1])
# number_color = "black"
# gaps_row = seq(4, 27, 4)
# cutree_col = 6
# scale = "column"

# Clear the working space so changes in the variables will take effect
rm(list = ls())

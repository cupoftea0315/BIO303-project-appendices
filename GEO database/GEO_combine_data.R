# Read data
rm(list=ls())
rawdata_combined_1 <- read.table("PRJNA615032 (GSE147507)/GSE147507_RawReadCounts_Human.tsv", header = TRUE, row.names = 1, fill = TRUE, sep="\t", na.strings=c("NA", ""))
rawdata_combined_2 <- read.table("PRJNA634167 (GSE150962)/GSE150962_merged_counts.txt", header = TRUE, row.names = 1, fill = TRUE, sep="\t", na.strings=c("NA", ""))

# Combine the data
rawdata_combined <- merge(rawdata_combined_1, rawdata_combined_2, by = "row.names", all = FALSE)

# Set row names
rownames(rawdata_combined) <- rawdata_combined$Row.names

# Remove the Row.names column
rawdata_combined$Row.names <- NULL

# Remove rows with NA values
rawdata_combined <- na.omit(rawdata_combined)

# Filter the dataset
selected_columns <- c(1:6, 84:89)
rawdata_combined <- rawdata_combined[, selected_columns]

# Define group variable
group <- c(rep("Normal", 3), 
           rep("SARS-CoV-2 infection", 6), 
           rep("Normal", 3))

# Determine batch based on column names
batch <- c(rep("GSE147507", 6), rep("GSE150962", 6))

# Save data as Rdata
save(rawdata_combined, group, batch, file = "GEO.Rdata")
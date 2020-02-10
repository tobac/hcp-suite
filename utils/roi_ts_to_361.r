#!/usr/bin/Rscript
# Reads a subject's time series files and combines them with Glasser parcels' time series
# Usage: $0 <input dir with all time series *txt files> <output_dir>
library(data.table)

# Assume we want to read all *.txt files as time series
args <- commandArgs(trailingOnly=TRUE)
input_dir <<- args[1]
output_dir <<- args[2]
data_dir <- gsub("HCP", "", Sys.getenv("HCP"))
parcels_dir <- "home/tobac/Projekte/diss/SCA21/glasser360"
file_list = list.files(path = input_dir, pattern = "*txt")
#file_list = readLines(args[1])

#dir.create(output_dir)

print("Reading ROI time series of all subjects...")
df_list_roi = lapply(file.path(input_dir, file_list), fread)

# Remove all subjects with incomplete data, i.e. less than 4800 time points
for(i in 1:length(df_list_roi)) { 
  if(length(df_list_roi[[i]]$V1) < 4800) { 
    df_list_roi[[i]] <- NA 
    file_list[[i]] <- NA
  } 
}

df_list_roi_clean = df_list_roi[!is.na(df_list_roi)] # Not needed anymore as we read txt files again later on
file_list_clean = file_list[!is.na(file_list)]
print(paste("Number of subjects with complete data:", length(file_list_clean)))

# Finally, combine ROI-time-series data frames with Glasser-parcellation-time-series data frames
#dir.create(output_dir)
options(scipen=10)
oldwd <- getwd()
for(i in 1:length(file_list_clean)) { # Yes, "for" loops are not the most elegant solution in R. They are easy to understand, though, and data integrity (combining identical IDs) is critical
  df_roi = fread(file.path(input_dir, file_list_clean[[i]])) # Seems superfluous as we already have this read in, but just to be sure we are combining corresponding files
  setwd(file.path(data_dir, parcels_dir))
  df_parcels = fread(basename(file_list_clean[[i]]))
  setwd(oldwd)
  df <- cbind(df_roi, df_parcels)
  output_file <- file.path(output_dir, basename(file_list_clean[[i]]))
  print(paste("Writing combined file", output_file, "|", i, "of", length(file_list_clean), "..."))
  write.table(df, file = output_file, row.names = FALSE, col.names = FALSE)
}

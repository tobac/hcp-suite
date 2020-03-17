ids <- scan("ids.file")
vars <- read.delim("vars.file")
vars <- as.character(read.table("vars.file", as.is = TRUE))
vars <- c(vars, "Subject")
# If you supply a string of space-separated variables instead of a file containing
# that string, conversion into a character vector of length n involves strsplit and thus
# uncommenting the following line
#vars <- strsplit(vars, split = " ")[[1]]

data1 <- read.csv("data1.csv")

if (file.exists("data2.csv")) {
  data2 <- read.csv("data2.csv")
  cat("Merging tables...\n")
  merged <- merge(data1, data2, by="Subject")
  write.csv(merged, file = "merged.csv")
} else {
  merged <- data1
}

if (file.exists("groups")) {
  groups <- scan("groups", what = "")
} else {
  # If no groups are specified, create one group spanning the entire sample
  groups <- paste("1:", ncol(merged)) 
}

cat("Picking variables and subjects...\n")
doi <- merged[merged$Subject %in% ids, ] # First get IDs of interest
doi <- doi[vars] # Then get columns of interest
if ("Gender" %in% names(doi)) {
  doi$Gender <- ifelse(doi$Gender == 'M',1,0) # We need to binarise gender
}

missing_values <- doi$Subject[which(rowSums(is.na(doi)) > 0)]
if (length(missing_values) > 0) {
  cat("\n")
  cat("The following subjects have missing values in at least one 
      of the specified variables.\n")
  print(missing_values)
  cat("Remove those subjects first\n")
  cat("\n")
  stop()
}

doi$Subject <- NULL # Remove Subject column

# Add group column(s)
for (idx in 1:length(groups)) { 
  colname <- paste("Z", idx, sep = "")
  doi[, as.character(colname)] <- 0 # Add group column and exclude everyone
  group <- strsplit(groups[idx], ":")[[1]][1]:strsplit(groups[idx], ":")[[1]][2]
  doi[group, as.character(colname)] <- 1 # Include members only
}

outfile = "design.mat"
cat("Writing design matrix...\n")
write(paste("/NumWaves", ncol(doi)), outfile)
write(paste("/NumPoints", nrow(doi)), outfile, append = TRUE)
ppheights <- sapply(doi, max) - sapply(doi, min)
write(paste("/PPHeights", paste(unlist(ppheights), collapse= " ")), outfile, append = TRUE)
write("", outfile, append = TRUE)
write("/Matrix", outfile, append = TRUE)
#options(scipen = -10) # We want scientific notation in our outfile
write.table(doi, file = outfile, row.names = FALSE, col.names = FALSE, append = TRUE)
options(scipen = 0)

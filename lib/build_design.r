data1 <- read.csv("data1.csv")

if (file.exists("data2.csv")) {
  data2 <- read.csv("data2.csv")
  cat("Merging tables...\n")
  merged <- merge(data1, data2, by="Subject")
  write.csv(merged, file = "merged.csv")
} else {
  merged <- data1
}

ids <- scan("ids.file")
vars <- read.delim("vars.file")
vars <- as.character(read.table("vars.file", as.is = TRUE))
vars <- c(vars, "Subject")
# If you supply a string of space-separated variables instead of a file containing
# that string, conversion into a character vector of length n involves strsplit and thus
# uncommenting the following line
#vars <- strsplit(vars, split = " ")[[1]]

cat("Picking variables and subjects...\n")
doi <- merged[merged$Subject %in% ids, ] # First get IDs of interest
doi <- doi[vars] # Then get columns of interest
if ("Gender" %in% names(doi)) {
  doi$Gender <- ifelse(doi$Gender == 'M',1,0) # We need to binarise gender
}

missing_values <- doi$Subject[which(rowSums(is.na(doi)) > 0)]
if (length(missing_values) > 0) {
  cat("\n")
  cat("The following subjects have missing values in of the specified variables.\n")
  print(missing_values)
  cat("Remove those subjects first\n")
  cat("\n")
  stop()
}
doi$Subject <- NULL # Remove Subject column
doi$Group <- 1 # Add group column

outfile = "design.mat"
cat("Writing design matrix...\n")
write(paste("/NumWaves", ncol(doi)), outfile)
write(paste("/NumPoints", nrow(doi)), outfile, append = TRUE)
write("", outfile, append = TRUE)
write("/Matrix", outfile, append = TRUE)
#options(scipen = -10) # We want scientific notation in our outfile
write.table(doi, file = outfile, row.names = FALSE, col.names = FALSE, append = TRUE)
options(scipen = 0)

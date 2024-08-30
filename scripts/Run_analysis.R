source("scripts/model_selection_by_var.R")
analyze_dependent_variable("Ganglia.Area.µm.2.")
analyze_dependent_variable("HuCD.Coverage")
analyze_dependent_variable("nNOS.Coverage")
analyze_dependent_variable("ChAT.Coverage")
analyze_dependent_variable("VIP.Coverage")
analyze_dependent_variable("SubP.Coverage")
analyze_dependent_variable("CGRP.Coverage")
analyze_dependent_variable("SNAP25.low.threshold.coverage")
analyze_dependent_variable("SNAP25.high.threshold.coverage")
analyze_dependent_variable("TH.Coverage")

library(openxlsx)

# Define the directory containing the CSV files
files <- list.files(path = "results/results_tables/", pattern = "*_model_results.csv", full.names = TRUE)
# Create a new workbook
wb <- createWorkbook()
# Loop through each file
for (file in files) {
  # Extract the sheet name from the filename
  sheet_name <- gsub("_model_results.csv$", "", basename(file))
  # Read the CSV file
  data <- read.csv(file)
  # Add a new sheet to the workbook
  addWorksheet(wb, sheet_name)
  # Write the data to the corresponding sheet
  writeData(wb, sheet = sheet_name, x = data)
}
# Save the workbook
saveWorkbook(wb, file = "results/results_tables/combined_results.xlsx", overwrite = TRUE)

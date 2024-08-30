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

#Excel
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

#Powerpoint
# Define the main directory containing the subdirectories with PNG files
main_dir <- "results/plots/"

# List all subdirectories in the main directory
subdirs <- list.dirs(main_dir, recursive = FALSE)

# Initialize a new PowerPoint object
ppt <- read_pptx()

# Loop through each subdirectory
for (subdir in subdirs) {
  # Get the cleaned dependent variable name from the subdir path
  cleaned_dependent_variable <- basename(subdir)
  
  # Add a title slide for the set of plots from this subdirectory
  ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = fpar(ftext(paste("Effect of predictor variables on", cleaned_dependent_variable),
                               fp_text(bold = TRUE, font.size = 24))), 
            location = ph_location_type(type = "title"))
  
  # List all PNG files in the current subdirectory
  png_files <- list.files(subdir, pattern = "\\.png$", full.names = TRUE)
  
  # Loop through each PNG file and add it as a slide
  for (png_file in png_files) {
    ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme") %>%
      ph_with(value = external_img(png_file), location = ph_location_fullsize())
  }
}

# Define the output PowerPoint file path in the main directory
output_pptx <- paste0(main_dir, "combined_plots.pptx")

# Save the PowerPoint presentation
print(ppt, target = output_pptx)

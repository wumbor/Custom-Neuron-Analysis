# This script analyses the data from optimise_threshold.ijm Fiji macro 
#Written by Victor Kumbol, June 2020
#Modified by Victor Kumbol, September 2024


##  Load required libraries and functions
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(ggrepel)


#list of treatment names
TreatmentLevels <- c("Null", "Control oligo", "Control oligo(10)","Alexa488", "miR-124-5p", 
                     "miR-124-5p(1)", "miR-124-5p(3)", "miR-124-5p(5)", "miR-124-5p(10)", "miR-124-5p(20)", 
                     "miR-9-5p", "miR-9-5p(1)", "miR-9-5p(3)", "miR-9-5p(5)", "miR-9-5p(10)", "miR-9-5p(20)", 
                     "miR-501-3p", "miR-501-3p(1)","miR-501-3p(3)", "miR-501-3p(5)", "miR-501-3p(10)", "miR-501-3p(20)",
                     "miR-92a-1-5p", "miR-92a-1-5p(1)", "miR-92a-1-5p(3)", "miR-92a-1-5p(5)", "miR-92a-1-5p(10)", 
                     "miR-92a-1-5p(20)", "let7b", "LOX", "R848", "TL8-506")

#General formatting for all graphs
plottheme <- theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) 


#FUNCTIONS
parse_standard_name <- function(input_string) {
  #function to parse common names of treatments to standard names
  
  codeNameVector <-  c("miR#4", "miR#5", "miR#13", "miR#19", "miR#27", "TL8", "miR92a", "miR92","miR124", "Alexa")
  standardNameVector <-  c("Control oligo", "miR-124-5p", "miR-9-5p", "miR-501-3p", "miR-92a-1-5p", "TL8-506", "miR-92a-1-5p", "miR-92a-1-5p", "miR-124-5p", "Alexa488")
  codeDoseVector <-  c("\\(L\\)", "\\(LM\\)", "\\(M\\)", "\\(H\\)",  "\\(XH\\)")
  standardDoseVector <- c("(1)", "(3)", "(5)", "(10)",  "(20)")
  
  output <- input_string
  i = 1
  for(i in 1:length(codeNameVector)){
    output <- str_replace(output, codeNameVector[i], standardNameVector[i])
  }
  
  i = 1
  for(i in 1:length(codeDoseVector)){
    output <- str_replace(output, codeDoseVector[i], standardDoseVector[i])
  }
  return(output)
}


determine_optimal_thresholds <- function(sequence_data, threshold_analysis_file) {
  
  ##return the minimum threshold determined for the control group
  thresholdHeaderLine <- (grep("Filename", read_lines(threshold_analysis_file, n_max = -1L))) - 1
  thresholdData <- read.csv(file = threshold_analysis_file, header = TRUE, sep = ",", skip = thresholdHeaderLine)
  
  thresholdData <- thresholdData %>%
    arrange(Filename) %>%
    mutate(Filename=trimws(Filename))
  
  sequence_data <- sequence_data %>%
    arrange(sequence_data) %>%
    mutate(sequence_data, Filename=trimws(str_remove(Filename, ".vsi")))
  
  thresholdData <- inner_join(sequence_data, thresholdData, by="Filename")  
  
  autoThreshold <- thresholdData %>%
    separate(Condition, sep = "_", into = c("Treatment", "Field")) %>%
    mutate(Field = as.numeric(Field)) %>%
    mutate(Treatment = parse_standard_name(trimws(Treatment))) %>%
    filter(str_detect(Treatment, "Null")) %>%
    summarise(NeuN.Threshold=min(NeuN.LowerThreshold), TUNEL.Threshold=min(TUNEL.LowerThreshold))

  return(autoThreshold)
}



#function to import microscopy sequence data
importCsvData <- function(filename, header) {
  header_line <- (grep(header, read_lines(filename, skip = 0, skip_empty_rows = FALSE, n_max = -1L, na = character()))) - 1 #Specify header line 
  df <- read.csv(file = filename, header = TRUE, sep = ",", skip = header_line)
  return(df)
}





analyze_threshold_data <- function(experiment_folder) {
  
  original_working_directory <- getwd()
  setwd(experiment_folder)

  experiment_id <- str_split(experiment_folder, "/", simplify = TRUE)
  experiment_id <- experiment_id[length(experiment_id)]
  
  #Specify input files
  microscopy_sequence_file <- paste0(experiment_id, "_Image_Capture.txt")
  sequence_data <- importCsvData(microscopy_sequence_file, "Condition")
  threshold_analysis_file <- paste0(experiment_id, "_ThresholdAnalysis.csv") 
  
  #Specify output files
  results_excel_file <- paste0(experiment_id, "_Optimised_Thresholds.xlsx")
  
  #determine optimal threshold  
  autoThreshold <- determine_optimal_thresholds(sequence_data, threshold_analysis_file)
  write.xlsx(autoThreshold, results_excel_file)
  
  setwd(original_working_directory)
}

####################################################################################################


#initialise workspace and allow user to define working directory
experiment_folder <- normalizePath(choose.dir(default = getwd(), caption = "Select Experiment Folder"), "/")
analyze_threshold_data(experiment_folder)




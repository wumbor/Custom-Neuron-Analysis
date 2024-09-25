# This script analyses the data from measure_neurite_properties_in_vivo Fiji macro 
# Written by Victor Kumbol, February 2023
# Modified by Victor Kumbol, September 2024

#load required libraries
library(tidyverse)
library(openxlsx)
library(ggpubr)
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





is_outlier <- function(x) {
  #function to test for outliers
  
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}




import_csv_data <- function(filename, header) {
  #function to import microscopy sequence data
  
  header_line <- (grep(header, read_lines(filename, skip = 0, skip_empty_rows = FALSE, n_max = -1L, na = character()))) - 1 #Specify header line 
  df <- read.csv(file = filename, header = TRUE, sep = ",", skip = header_line)
  return(df)
}


import_neurite_data <- function(microscopyDataFile, imageSequence) {
  #function to import neurite analysis data
  
  #read in data from files
  imageData <- read.csv(file = microscopyDataFile, header = TRUE, sep = ",")
  
  imageSequence <- imageSequence %>% #process image sequence data
    rename(Original_Filename = Filename) %>%
    arrange(Original_Filename)
  
  imageData <- imageData %>%  #process image data
    rename(Processed_Filename = Filename) %>%
    arrange(Processed_Filename) %>%
    mutate(Degeneration.Index = FragmentedNeuriteArea/TotalNeuriteArea) %>%
    mutate(NeuriteLengthPerCell=TotalNeuriteLength/NucleiCount) %>%
    mutate(NeuriteAttachmentPointsPerCell=NeuriteAttachmentPoints/NucleiCount)    
  
  importedData <- bind_cols(imageSequence, imageData) %>% #merge sequence and data tables 
    separate(Condition, sep = "_", into = c("Treatment", "Field")) %>%
    mutate(Field = as.numeric(Field)) %>%
    mutate(Treatment = parse_standard_name(trimws(Treatment)))
  importedData$Treatment <- factor(importedData$Treatment, levels = TreatmentLevels)
  return(importedData)
}


process_neurite_data <- function(DataToProcess, resultsExcelFile) {
  
  #function to process neurite data
  
  if (file.exists(resultsExcelFile)){   #overwrite any existing results file
    file.remove(resultsExcelFile)
  } 
  
  DataToProcess <- DataToProcess %>% #identify suspected outliers to aid visual inspection 
    filter(!is.na(Degeneration.Index)) %>%
    filter(NucleiCount>0) %>% #exclude fields without any cells
    mutate(Suspected.Outliers.Degeneration.Index=Processed_Filename) %>% 
    mutate(Suspected.Outliers.NeuriteLengthPerCell=Processed_Filename) %>% 
    group_by(Treatment) %>% 
    mutate(is_outlier.degeneration.index=ifelse(is_outlier(Degeneration.Index), Degeneration.Index, as.numeric(NA))) %>%
    mutate(is_outlier.neurite.length=ifelse(is_outlier(NeuriteLengthPerCell), NeuriteLengthPerCell, as.numeric(NA)))
  DataToProcess$Suspected.Outliers.Degeneration.Index[which(is.na(DataToProcess$is_outlier.degeneration.index))] <- as.numeric(NA)
  DataToProcess$Suspected.Outliers.NeuriteLengthPerCell[which(is.na(DataToProcess$is_outlier.neurite.length))] <- as.numeric(NA)
  DataToProcess <- DataToProcess %>%
    select(-c(is_outlier.degeneration.index, is_outlier.neurite.length))
  
  summaryReport <- DataToProcess %>% #Generate summary results
    group_by(Treatment) %>%
    summarise(Mean.Degeneration.Index = mean(Degeneration.Index), Mean.NeuriteLengthPerCell = mean(NeuriteLengthPerCell))
  
  #save the results to an excel spreadsheet
  data_list <- list("Raw Data"=DataToProcess, "Neurite Analysis Summary"=summaryReport)
  write.xlsx(data_list, resultsExcelFile, col.names = TRUE, row.names = FALSE)
  
  return(DataToProcess)
}



plot_neurite_data <- function(DataToPlot, resultsGraphFile, experiment_id) {
  
  #function to plot neurite data
  
  scaledWidth = 550*length(unique(DataToPlot$Treatment)) #scale width to accommodate treatments
  
  p1 <- ggplot(DataToPlot, aes(y=NeuriteLengthPerCell, x=Treatment, fill = factor(Mouse.ID))) + 
    geom_boxplot() +  stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") + 
    geom_text_repel(aes(label = Suspected.Outliers.Degeneration.Index), na.rm = TRUE, show.legend = FALSE) + 
    theme_pubr() + plottheme
  
  p2 <- ggplot(DataToPlot, aes(y=Degeneration.Index, x=Treatment, fill = factor(Mouse.ID))) + 
    geom_boxplot() +  stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") + 
    geom_text_repel(aes(label = Suspected.Outliers.Degeneration.Index), na.rm = TRUE, show.legend = FALSE) + 
    theme_pubr() + plottheme
  
  p <- ggarrange(plotlist = list(p1, p2), common.legend = T)
  
  Graphs <- annotate_figure(p, top = text_grob(experiment_id, face = "bold", size = 14))
  ggsave(resultsGraphFile, Graphs, width = scaledWidth*2, height = 3750, dpi = 300, units = "px",  compression="lzw")
  
}




analyse_neurite_properties <- function(experiment_folder) {
  #function to execute the complete analysis
  
  original_working_directory <- getwd()
  setwd(experiment_folder)
  
  #define key global variables
  experiment_id <- str_split(experiment_folder, "/", simplify = TRUE)
  experiment_id <- experiment_id[length(experiment_id)]

  #Specify input files
  microscopy_sequence_file <- paste(experiment_id, "_Image_Capture.txt", sep = "")
  sequence_data <- import_csv_data(microscopy_sequence_file, "Condition")
  
  microscopyDataFile <- paste(experiment_id, "_", neurite_marker, "_Neurite_Analysis.csv", sep = "")
  resultsExcelFile <- paste(experiment_id, "_",neurite_marker, "_Neurite_Degeneration.xlsx", sep = "")
  resultsGraphFile <- paste(experiment_id, "_",neurite_marker, "_Neurite_Degeneration.tiff", sep = "")
  
  DataToProcess <- import_neurite_data(microscopyDataFile, sequence_data) #import data
  DataToPlot <- process_neurite_data(DataToProcess, resultsExcelFile) #process data
  plot_neurite_data(DataToPlot, resultsGraphFile, experiment_id) #plot graphs
  
  setwd(original_working_directory)
  
}

####################################################################################################


experiment_folder <- normalizePath(choose.dir(default = getwd(), caption = "Select Experiment Folder"), "/")
analyse_neurite_properties(experiment_folder)





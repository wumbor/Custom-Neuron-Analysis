# This script analyses the data output of the miRKinetics macro to generate miRNA uptake kinetics graphs
# Written by Victor Kumbol, August 2021

#load required libraries
library(tidyverse)
library(ggpubr)
library(ggsci)
library(openxlsx)
library(ggplot2)


#define custom functions
extract_image_parameter <- function(parameter_name, image_parameters) { 
  parameter <- image_parameters$Value[str_detect(image_parameters$Parameter, parameter_name)]
  parameter <- trimws(parameter)
  parameter <- str_split(parameter, " ", simplify = TRUE)
  parameter <- parameter[1]
  return(parameter)
}


import_kinetics_data  <- function(image_data_folder) {
  original_working_directory <- getwd()
  setwd(image_data_folder)
  files_to_process <- list.files(image_data_folder, pattern = "*.csv") #retrieve list of results files
  imported_data <- list()
  
  for (i in 1:length(files_to_process)) { #iterate over the files in the folder and import each
    raw_data <- read.csv(files_to_process[i])
    raw_data <- raw_data %>%
      select(Channel, Area, Mean, RawIntDen) %>%
      separate(Channel, into = c("Channel", "Slice"), sep = "-") %>%
      mutate(Channel = trimws(Channel))
    #append results to the imported_data list
    imported_data <- append(imported_data, list(raw_data))
  }
  names(imported_data) <- files_to_process
  
  #merge results into one
  imported_data <- bind_rows(imported_data, .id="Filename")
  imported_data <- imported_data %>%
    rename(Parameter = Channel) %>%
    mutate(Parameter = str_replace(Parameter, "red inverse", "redInverse")) %>%
    separate(Filename, into = c("Filename", "Channel.FrameID"), sep = "_miRKinetics_" ) %>%
    separate(Channel.FrameID, into = c("Channel", "FrameID"), sep = "-") %>%
    separate(FrameID, into = c("FrameID", NA), sep = "\\.") %>%
    mutate(FrameID = as.numeric(FrameID), Slice = as.numeric(Slice))
  
  setwd(original_working_directory)
  return(imported_data)
}



process_kinetics_data <- function(combined_data, time_offset, frame_interval, slice_selection, resultsExcelFile, parameter_labels, channel_labels) { #This function processes the imported data table and exports them
  
  if (file.exists(resultsExcelFile)){   #overwrite any existing results file
    file.remove(resultsExcelFile)
  } 
  
  combined_data_per_slice <- combined_data %>%  
    filter(Slice %in% slice_selection) %>%
    pivot_wider(names_from = Parameter, values_from = c(Mean, Area, RawIntDen)) %>% 
    mutate(Mean_prop = Mean_red/(Mean_red + Mean_redInverse), RawIntDen_prop = RawIntDen_red/(RawIntDen_red + RawIntDen_redInverse), Area_prop = Area_red/(Area_red + Area_redInverse)) %>% #compute other parameters
    filter(!is.na(Area_redInverse)) %>% #filter out slices where the whole area was thresholded
    select(-Area_redInverse, -RawIntDen_redInverse, -Mean_redInverse) %>%
    pivot_longer(Mean_red:Area_prop, names_to = "Parameter", values_to = "Value") %>%
    mutate(RecordingTime.h = ((FrameID-1)*frame_interval)/60) %>%
    mutate(AbsoluteTime.h = RecordingTime.h + (time_offset/60)) %>%
    rename(Mean = Value)
   
  combined_data_per_slice <- left_join(combined_data_per_slice, parameter_labels, by = "Parameter")
  combined_data_per_slice <- left_join(combined_data_per_slice, channel_labels, by = "Channel")
  
  
  combined_data_per_frame <- combined_data_per_slice %>% #summarise data for each frame
    select(-c(Parameter, Channel, RecordingTime.h)) %>%
    group_by(Filename, AbsoluteTime.h, ChannelLabel, Region) %>%
    summarise(Mean = mean(Mean))
  combined_data_per_frame <- combined_data_per_frame %>% #calculate normalised values
    group_by(Filename, ChannelLabel, Region) %>%
    mutate(NormalisedMean = (Mean/first(Mean))) %>%
    ungroup()
    
  write.xlsx(as.data.frame(combined_data_per_frame), file = resultsExcelFile, sheetName = "combinedDataPerFrame", col.names = TRUE, row.names = FALSE)
  
  return(combined_data_per_slice)
}







plot_mir_kinetics_graphs <- function(plot_data, duration_for_plot, errorType, treatment, resultsGraphFile) {

  plot_data <- plot_data %>%
    filter(AbsoluteTime.h <= duration_for_plot) #plot only 4h results
  
  #Plot graph of Endosomal Mean
  p1 <- plot_data %>%
    filter(Region=="Endosomal Signal Mean") %>%
    ggline(x= "AbsoluteTime.h", y = "Mean", ylab = "Endosomal Signal - Mean", color ="ChannelLabel", palette =  c("#0000FF", "#00FF00", "#FF0000"), numeric.x.axis = TRUE, xticks.by = 1, title =treatment, add = errorType)
  
  ggexport(p1, filename = resultsGraphFile, width = 1500, height = 750, res=200, compression="lzw")

}




analyze_kinetics_experiment <- function(experiment_folder) {
  
  original_working_directory <- getwd()
  setwd(experiment_folder)
  
  #define key global variables
  experiment_id <- str_split(experiment_folder, "/", simplify = TRUE)
  experiment_id <- experiment_id[length(experiment_id)]
  
  #Specify input files
  image_parameters_file <- paste0(experiment_id, "_Image_Capture.txt")
  image_data_folder <- paste0(experiment_folder, "/", experiment_id, "_Results")
  header_line <- (grep("Parameter", read_lines(image_parameters_file, skip = 0, skip_empty_rows = FALSE, n_max = -1L, na = character()))) - 1
  
  #Read analysis parameters from file
  imageParameters <- read.csv(file = image_parameters_file, header = TRUE, sep = ":", skip = header_line)
  treatment <- extract_image_parameter("Treatment", imageParameters)
  frame_interval <- as.numeric(extract_image_parameter("Frame Interval", imageParameters))
  time_offset <- as.numeric(extract_image_parameter("Time Offset", imageParameters))
  slice_selection_range <- extract_image_parameter("Slice Selection", imageParameters)
  slice_selection_range <- str_split(slice_selection_range, "-", simplify = TRUE)
  slice_selection <- seq(slice_selection_range[1], slice_selection_range[2],by=1)
  
  #define designations for parameters
  channel_labels<- tibble(Channel = c("red", "green", "blue"), ChannelLabel = c("phRodored", "miRNA", "DAPI"))
  parameter_labels <- tibble(Parameter = c("Mean_red", "Mean_prop", "Area_red", "Area_prop", "RawIntDen_red", "RawIntDen_prop"), Region = c("Endosomal Signal Mean","Endosomal Mean Proportion","Endosomal Area",  "Endosomal Area Proportion", "Endosomal Signal RawIntDen", "Endosomal RawIntDen Proportion"))
  resultsExcelFile <- paste(treatment, "_Kinetics.xlsx", sep = "")
  resultsGraphFile <- paste(treatment, "_Kinetics.tiff", sep = "")
  
  
  #run data analysis
  kinetics_data  <- import_kinetics_data (image_data_folder)
  plot_data <- process_kinetics_data(kinetics_data, time_offset, frame_interval, slice_selection, resultsExcelFile, parameter_labels, channel_labels)
  
  #plot graphs
  plot_mir_kinetics_graphs(plot_data, duration_for_plot = 4, errorType = "mean_sd", treatment, resultsGraphFile)
  setwd(original_working_directory)
}




####################################################################################################


#initialise workspace and allow user to define working directory
experiment_folder <- normalizePath(choose.dir(default = getwd(), caption = "Select Experiment Folder"), "/")
analyze_kinetics_experiment(experiment_folder)


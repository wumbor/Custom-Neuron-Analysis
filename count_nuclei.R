# This script analyses the data from the count_nuclei.ijm Fiji macro 
#Written by Victor Kumbol, June 2020
#Modified by  Victor Kumbol, September 2024 

##  Load required libraries and functions
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(ggrepel)


#list of treatment names
TreatmentLevels <- c("Null", "Control oligo","Control oligo(5)", "Control oligo(10)","Alexa488", "miR-124-5p", 
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





####################################################################################################


#function to import microscopy sequence data
importCsvData <- function(filename, header) {
  header_line <- (grep(header, read_lines(filename, skip = 0, skip_empty_rows = FALSE, n_max = -1L, na = character()))) - 1 #Specify header line 
  df <- read.csv(file = filename, header = TRUE, sep = ",", skip = header_line)
  return(df)
}


normalize <- function(value, maxi) {
  normalized_value = (value/maxi)*100
  return(normalized_value)
}

#function to test for outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}


import_count_data <- function(microscopy_data_file, sequence_data) {
  #read in data from files
  imageData <- read.csv(file = microscopy_data_file, colClasses = c("character", rep("numeric", 4)), header = TRUE, sep = ",")
  imageData <- imageData %>%  #process image count data
    separate(Slice, into = c("Channel", "Filename"), sep = "-") %>%
    select(Channel, Filename, Count)  %>%
    mutate(Filename=trimws(str_remove(Filename, ".tif")))
  
  sequence_data <- sequence_data %>%
    arrange(sequence_data) %>%
    mutate(sequence_data, Filename=trimws(str_remove(Filename, ".vsi")))
  
  imageData <- inner_join(sequence_data, imageData, by="Filename")  
  
  imageData <- imageData %>% #process image sequence data
    separate(Condition, sep = "_", into = c("Treatment", "Field")) %>%
    mutate(Field = as.numeric(Field)) %>%
    mutate(Treatment = parse_standard_name(trimws(Treatment)))
  imageData$Treatment <- factor(imageData$Treatment, levels = TreatmentLevels)
  return(imageData)
}






process_count_data <- function(DataToProcess, resultsExcelFile, NeuNChannel, TUNELChannel) {
  if (file.exists(resultsExcelFile)){   #overwrite any existing results file
    file.remove(resultsExcelFile)
  } 

  DataToProcess <- DataToProcess %>%
    mutate(Channel=str_replace(Channel, NeuNChannel, "NeuN.Count")) %>%
    mutate(Channel=str_replace(Channel, TUNELChannel, "TUNEL.Count")) %>%
    pivot_wider(names_from = Channel, values_from = Count)
    
  DataToProcess <- DataToProcess %>%  
    mutate(Suspected.Outliers.NeuN.Count=Filename) %>% 
    mutate(Suspected.Outliers.TUNEL.Count=Filename) %>% 
    group_by(Treatment) %>% 
    mutate(is_outlier.NeuN.Count=ifelse(is_outlier(NeuN.Count), NeuN.Count, as.numeric(NA))) %>%
    mutate(is_outlier.TUNEL.Count=ifelse(is_outlier(TUNEL.Count), TUNEL.Count, as.numeric(NA)))
  
  DataToProcess$Suspected.Outliers.NeuN.Count[which(is.na(DataToProcess$is_outlier.NeuN.Count))] <- as.numeric(NA)
  DataToProcess$Suspected.Outliers.TUNEL.Count[which(is.na(DataToProcess$is_outlier.TUNEL.Count))] <- as.numeric(NA)
  DataToProcess <- DataToProcess %>%
    select(-c(is_outlier.NeuN.Count, is_outlier.TUNEL.Count)) %>%
    mutate(TUNEL_NeuN_Ratio = TUNEL.Count/NeuN.Count)
  
  summaryReport <- DataToProcess %>%
    group_by(Treatment) %>%
    summarise(NeuNCount=mean(NeuN.Count), TUNELCount=mean(TUNEL.Count)) %>%
    mutate(TUNEL_NeuN_Ratio = TUNELCount/NeuNCount)
  
  nullGroupStats <- summaryReport %>% #calculate mean of null group
    filter(Treatment=="Null") 
  
  summaryReport <- summaryReport %>% #Normalize all counts to null group stats
    mutate(Normalized.NeuNCount = normalize(NeuNCount, nullGroupStats$NeuNCount))  %>%
    mutate(Normalized.TUNEL_NeuN_Ratio = normalize(TUNEL_NeuN_Ratio, nullGroupStats$TUNEL_NeuN_Ratio)/100)
  
  #save the results to an excel spreadsheet
  data_list <- list("Raw Data"=DataToProcess, "Nuclei Count Summary"=summaryReport)
  write.xlsx(data_list, resultsExcelFile, col.names = TRUE, row.names = FALSE)
  
  return(DataToProcess)
}


plotCountGraphs <- function(DataToPlot, resultsGraphFile) {
   scaledWidth = 450*length(unique(DataToPlot$Treatment)) #scale width to accommodate treatments
  
  p1 <- ggplot(DataToPlot, aes(y=NeuN.Count, x=Treatment)) + 
    geom_boxplot() +  stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") + 
    geom_text_repel(aes(label =Suspected.Outliers.NeuN.Count), na.rm = TRUE, show.legend = FALSE) + 
    theme_pubr() + plottheme
  
  p2 <- ggplot(DataToPlot, aes(y=TUNEL.Count, x=Treatment)) + 
    geom_boxplot() +  stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") + 
    geom_text_repel(aes(label = Suspected.Outliers.TUNEL.Count), na.rm = TRUE, show.legend = FALSE) + 
    theme_pubr() + plottheme
  
  p <- ggarrange(plotlist = list(p1, p2), common.legend = T)
  Graphs <- annotate_figure(p, top = text_grob(experiment_id, face = "bold", size = 14))
  ggsave(resultsGraphFile, Graphs, width = scaledWidth*2, height = 3750, dpi = 300, units = "px",  compression="lzw")
}




analyze_nuclei_counts <- function(experiment_folder) {
  
  original_working_directory <- getwd()
  setwd(experiment_folder)

  experiment_id <- str_split(experiment_folder, "/", simplify = TRUE)
  experiment_id <- experiment_id[length(experiment_id)]
  
  #Specify input files
  microscopy_sequence_file <- paste0(experiment_id, "_Image_Capture.txt")
  sequence_data <- importCsvData(microscopy_sequence_file, "Condition")
  microscopy_data_file <- paste0(experiment_id, "_Nuclei_Count.csv")

  
  #Specify output files
  resultsExcelFile <- paste0(experiment_id, "_Nuclei_Count.xlsx")
  resultsGraphFile <- paste(experiment_id, "_Nuclei_Count.tiff")
  
  #read nuclei count parameters from file  
  DataToProcess <- import_count_data(microscopy_data_file, sequence_data)
  DataToPlot <- process_count_data(DataToProcess, resultsExcelFile, NeuNChannel, TUNELChannel)
  plotCountGraphs(DataToPlot, resultsGraphFile) #plot graphs
  
  setwd(original_working_directory)
}

####################################################################################################

#initialise workspace and allow user to define working directory
NeuNChannel <- "C1"
TUNELChannel <- "C2"
experiment_folder <- normalizePath(choose.dir(default = getwd(), caption = "Select Experiment Folder"), "/")
analyze_nuclei_counts(experiment_folder)




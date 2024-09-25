//This macro analyzes specified channels from a batch of images and returns a list of thresholds after running an autothreshold algorithm
//Written by Victor Kumbol, June 2020
//Modified by Victor Kumbol, September 2024


/*
 	***********************

	Set IJ parameters

	***********************
*/

setBatchMode(true);
run("Clear Results");

/*
 	***********************

	Set experiment variables

	***********************
*/
var workingDir = "C:/User/Desktop/NeuriteImages/ExperimentID"; 			//set experiment folder
var NeuNChannel = "C1-";												//channel for NeuN nuclei
var TUNELChannel = "C2-";												//channel for TUNEL nuclei

/*
 	***********************

		Functions

	***********************
*/

function getThresholdData(windowtitle, filename, thresholdAnalysisFile) { 
	//Extract the appropriate channel from image
	run("Split Channels");
	
	//Get NeuN lower threshold
	selectWindow(NeuNChannel + windowtitle);
	setAutoThreshold("Otsu dark");
	getThreshold(lowerThreshold, upperThreshold);
	NeuNThreshold=lowerThreshold;
	resetThreshold;
	
	//Get TUNEL lower threshold
	selectWindow(TUNELChannel + windowtitle);
	setAutoThreshold("Otsu dark");
	getThreshold(lowerThreshold, upperThreshold);
	TUNELThreshold=lowerThreshold;
	resetThreshold;	
	
	//write data to file
	File.append(filename + ", " + NeuNThreshold + ", " + TUNELThreshold, thresholdAnalysisFile);
	close("*");
}


function runThresholdAnalysis(workingDir) { 
// this function executes the analysis by calling other sub-functions
	experimentId = split(workingDir, "/"); 
	experimentId = experimentId[(lengthOf(experimentId)-1)]; //Extract the experimentId from the directory path
	sourceImagesDir = workingDir + "TIFFs/"; //define source images directory
	fileList = getFileList(sourceImagesDir);

	thresholdAnalysisFile = workingDir + experimentId + "_ThresholdAnalysis.csv";
	//Prepare the ThresholdAnalysis file 
	garbage = File.delete(thresholdAnalysisFile);	////delete any existing Threshold analysis file
	f = File.open(thresholdAnalysisFile); 
	print(f, "Filename, NeuN.LowerThreshold, TUNEL.LowerThreshold");	
	File.close(f);
	
	
	//Get threshold data for all images in the selected folder
	print("\n\nFinding the optimum threshold for nuclei count...");
	for (i = 0; i < fileList.length; i++){
		if (endsWith(fileList[i], ".tif")) { //process only tiff images
			filePath = sourceImagesDir + fileList[i];
			open(filePath);
			if (is("composite")) { //process only composite images
				getThresholdData(getTitle(), File.nameWithoutExtension, thresholdAnalysisFile);
			}
		}
	    close("*");  
	    showProgress(i, fileList.length);	
	}
	print("Threshold analysis complete");

}


/*
 	***********************

	Execute code to analyze images and save results

	***********************
*/

runThresholdAnalysis(workingDir);








//This macro counts the number of nuclei in images and saves the processed images, if requested by the user
//Written by Victor Kumbol, June 2020
//Modified by Victor Kumbol, September 2024


/*
 	***********************

	Set IJ parameters

	***********************
*/

setBatchMode(true);
run("Clear Results");
run("Set Measurements...", "area mean modal median limit display redirect=None decimal=3");



/*
 	***********************

	Set experiment variables

	***********************
*/

//Retrieve macro arguments
var workingDir = "C:/User/Desktop/NeuriteImages/ExperimentID"; 				//set experiment folder
var NeuNChannel = "C1-";													//channel for NeuN nuclei
var TUNELChannel = "C2-";													//channel for TUNEL nuclei
var NeuNThreshold = 485;													//optimised threshold for NeuN nuclei
var TUNELThreshold = 245;													//optimised threshold for TUNEL nuclei
var saveImages = true;														//should processed images be saved?


/*
 	***********************

		Functions

	***********************
*/

function countNeuNNuclei(resultsDir, fileName, NeuNChannelWindow) { 
	
	//Duplicate NeuN channel from image
	selectWindow(NeuNChannelWindow);
	run("Duplicate...", "title=original_neun_image");

	//Process image
	selectWindow(NeuNChannelWindow);
	run("Median...", "radius=5");
	run("Subtract Background...", "rolling=50");
	setThreshold(NeuNThreshold, 65536);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Fill Holes");
	run("Convert to Mask");
	run("Watershed");
	
	//Analyze image for nuclei
	run("Analyze Particles...", "size=20-Infinity show=Outlines exclude summarize");

	if (saveImages) {
		//Save masked image with counted cells
		selectWindow("original_neun_image");
		run("Grays");
		run("Apply LUT");
		run("8-bit");
		
		selectWindow("Drawing of " + NeuNChannelWindow);
		rename("neun_count_mask_image");
		run("Magenta");
		run("Invert LUT");
		run("Merge Channels...", "c4=&original_neun_image c6=&neun_count_mask_image");
		saveAs("tiff", resultsDir + fileName + "_nuclei_count_mask");
	}

}



function countTUNELNuclei(resultsDir, fileName, TUNELChannelWindow) { 
	
	//Duplicate TUNEL channel from image
	selectWindow(TUNELChannelWindow);;
	run("Duplicate...", "title=original_tunel_image");
	
	//Process image
	selectWindow(TUNELChannelWindow);
	run("Median...", "radius=5");
	run("Subtract Background...", "rolling=50");
	setThreshold(TUNELThreshold, 65536);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Fill Holes");
	run("Convert to Mask");
	run("Watershed");
	
	//Analyze image for nuclei
	run("Analyze Particles...", "size=10-Infinity show=Outlines exclude summarize");

	if (saveImages) {
		//Save masked image with counted cells
				
		//Adjust contrast of TUNL count images
		selectWindow("original_tunel_image");
		run("Grays");
		run("Apply LUT");
		run("Enhance Contrast", "saturated=0.35");
		run("Apply LUT");
		setOption("ScaleConversions", true);
		run("8-bit");

		selectWindow("Drawing of " + TUNELChannelWindow);
		rename("tunel_count_mask_image");
		run("Magenta");
		run("Invert LUT");
		run("Merge Channels...", "c4=&original_tunel_image c6=&tunel_count_mask_image");
		saveAs("tiff", resultsDir + fileName + "_tunel_count_mask");
	}

}


function deleteFolder(folderpath) { 
	//first delete all files in the folder
	deleteFileList = getFileList(folderpath);
		for (i = 0; i < deleteFileList.length; i++){
			deletefilePath = folderpath +deleteFileList[i];
			garbage = File.delete(deletefilePath);
			}	

	//then delete the folder itself
	deletecompleted = File.delete(folderpath);
	if (deletecompleted) {
		return true;
	} else {
		return false;
	}
}



function analyzeNucleiCounts(workingDir) { 
	// this function executes the analysis by calling other sub-functions
	
	//Define key variables
	experimentId = split(workingDir, "/"); 
	experimentId = experimentId[(lengthOf(experimentId)-1)]; //Extract the experimentId from the directory path
	sourceImagesDir = workingDir + "TIFFs/"; //define source images directory
	fileList = getFileList(sourceImagesDir);
	
	resultsDir = workingDir + "Nuclei Count/"; //define output images directory
	resultsFile = workingDir + experimentId + "_Nuclei_Count.csv";
			

	//Create a results subfolder
	if (!File.exists(resultsDir)) {
		File.makeDirectory(resultsDir);
		print("");
		print("\nResults folder created: ");
		print(resultsDir);
		print("");
	} else {
		if (deleteFolder(resultsDir)) {
			print("\nExisting Results Folder Deleted");
			File.makeDirectory(resultsDir);
			print("\nNew Results folder created: ");
		}
	}
		
	
	//Run the nuclei count function on all images in the selected folder
	print("\n\nCounting Nuclei...");
	for (i = 0; i < fileList.length; i++){
	//for (i = 0; i < 3; i++){
		if (!(endsWith(fileList[i], ".tif"))) {
			exit("Error: Please provide TIFF images for processing"); // return an error and exit the macro if the image is not a TIFF
		}
		filePath = sourceImagesDir + fileList[i];
		open(filePath);
		if (!(is("composite"))) {
			exit("Error: Please provide composite images for processing"); // return an error and exit the macro if the image is not multichannel
		}
		
		imageWindow = getTitle();
		fileName = File.nameWithoutExtension;
		run("Split Channels");
		NeuNChannelWindow = NeuNChannel+imageWindow;
		TUNELChannelWindow = TUNELChannel+imageWindow;	
		
		countNeuNNuclei(resultsDir, fileName, NeuNChannelWindow);
		countTUNELNuclei(resultsDir, fileName, TUNELChannelWindow);	
		close("*");	
		showProgress(i, fileList.length);	
	}
		
	
	selectWindow("Summary"); 
	saveAs("Results", resultsFile);
	run("Close All");
	
	print("Nuclei Count Complete");
	print("Total Files Processed: " + fileList.length);
	
}



/*
 	***********************

	Execute code to analyze images and save results

	***********************
*/

analyzeNucleiCounts(workingDir);
















//This macro measures neurite length and additional parameters from fluorescence images of mouse cortical brain sections. The sections should be stained for at least one neurite marker and one nuclear marker.
//Written by Victor Kumbol, August 2022
//Modified by Victor Kumbol, September 2024


/*
 	***********************

	Set IJ parameters

	***********************
*/
setBatchMode(false);
run("Clear Results");
run("Set Measurements...", "area limit display redirect=None decimal=3");
run("Conversions...", "scale calibrate");

/*
 	***********************

	Set experiment variables

	***********************
*/

var workingDir ="C:/User/Desktop/NeuriteImages/ExperimentID/";				//set experiment folder
var fixThreshold = false;													//should user-defined thresholds be applied?
var somaChannel = "C3-";													//channel for nuclei
var neuriteChannel = "C2-";													//channel for neurites
var somaThreshold = 0;														//user-defined threshold for nuclei
var neuriteThreshold = 0; 													//user-defined threshold for neurites
var neuriteMarker = "MAP2"; 												//name of neurite marker 
var saveNucleiMask = false;													//should nuclei mask image be saved?
var saveNeuriteMask = true;													//should neurite mask image be saved?
var saveFragmentationMask = true;											//should neurite fragments mask image be saved?

var pixel_size = 0;
var pixel_unit = "";



/*
 	***********************

		Functions

	***********************
*/

function prepareSomaImage(imageWindow) { 
// this function segments nuclei from the image

	//process image to enhance features for detection
	selectWindow(imageWindow);
	run("Subtract Background...", "rolling=50");
	run("Enhance Local Contrast (CLAHE)", "blocksize=100 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
	run("Gaussian Blur...", "sigma=3"); 
	
	//segment nuclei from image
	if (fixThreshold) {
		setThreshold(threshold, 65535);
	} else {
		setAutoThreshold("MaxEntropy dark");
	}
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Watershed");
	run("Fill Holes");
	
	//retain only nuclei > 20 square microns
	run("Analyze Particles...", "size=20-Infinity show=Masks");
	nucleiMaskImage = getTitle();
	selectWindow(nucleiMaskImage);
	run("Invert LUT");	
	rename("nuclei_mask");

}




function prepareNeuriteImage(imageWindow) { 
// this function segments neurites from the image

	getImageCalibration(imageWindow);
	selectWindow(imageWindow);
	run("Frangi Vesselness", "input=[imageWindow] dogauss=true spacingstring=[1, 1] scalestring=[2, 5]");
	run("Make Subset...", "slices=1");
	setAutoThreshold("Huang dark no-reset");
	run("Convert to Mask");
	
	//remove speckles
	run("Analyze Particles...", "size=0.20-Infinity show=Masks");
	run("Invert LUT");
	rename("neurite_mask");
	
	//recalibrate neurite mask image
	calibrateImage("neurite_mask");
	
	imageCalculator("Subtract create", "neurite_mask", "nuclei_mask");
	selectImage("Result of neurite_mask");
	rename("neurite without nuclei mask");
	
	//prepare neuronal network mask
	imageCalculator("Add create", "nuclei_mask", "neurite without nuclei mask");
	selectImage("Result of nuclei_mask");
	rename("neuronal network mask");
	
}





function measureNeuriteProperties(imageWindow, resultsFile, resultsDir) { 
// this function calculates various neurite properties (total length, number and area) given the image of a soma and dendrite

	//Image reconstruction and analysis
	selectImage("neurite_mask");
	run("Duplicate...", "title=[neurite_mask-skeleton]");
	run("Skeletonize");
	selectImage("nuclei_mask");
	run("Duplicate...", "title=[nuclei_mask-attachments]");
	run("Options...", "iterations=1 count=3 black do=Dilate");

	
	//Get neurite attachment points data
	imageCalculator("AND create", "nuclei_mask-attachments","neurite without nuclei mask");
	selectImage("Result of nuclei_mask-attachments");
	run("Find Maxima...", "prominence=0 strict output=Count");
	NeuriteAttachmentPoints = getResult("Count");
	run("Clear Results");
	
	
	//Get neurite attachment point mask
	selectImage("Result of nuclei_mask-attachments");
	run("Find Maxima...", "prominence=0 strict output=[Point Selection]");
	run("Flatten");
	selectImage("Result of nuclei_mask-attachments-1");
	rename("neurite attachment points mask");
	selectImage("neurite attachment points mask");
	run("8-bit");
	
	
	//Get neurite length data
	selectImage("neurite_mask-skeleton");
	run("Summarize Skeleton");
	selectWindow("Skeleton Stats");
	TotalNeuriteLength = getResult("Total length");
	MaxBranchLength = getResult("Max branch length");
	MeanBranchLength = getResult("Mean branch length");
	TotalNeuriteEndPoints = getResult("# End-points");
	TotalBranches = getResult("# Branches");
	TotalTrees = getResult("# Trees");
	run("Clear Results");
	
	
	//Get nuclei count
	selectImage("nuclei_mask");
	run("Analyze Particles...", "size=20-Infinity show=Outlines summarize");
	NucleiCount = getResult("Count");
	NucleiTotalArea = getResult("Total Area");
	NucleiAverageArea = getResult("Average Size");
	run("Clear Results");

	
	//Get total neuron area 
	selectImage("neuronal network mask");
	run("Create Selection");
	run("Measure");
	TotalNeuronArea = getResult("Area");
	run("Clear Results");
	close("Results");
	close("Summary");
	
	
	//Measure total neurite area
	selectImage("neurite_mask");
	run("Analyze Particles...", "size=0.0202-Infinity circularity=0-1.00 summarize");


	//Measure fragments area
	selectImage("neurite_mask");
	run("Analyze Particles...", "size=0.0202-0.255 circularity=0.2-1.00 show=[Bare Outlines] summarize");
	rename("fragment_mask_image");
	
	//Get total neurite and fragments area
	selectWindow("Summary");
	IJ.renameResults("Summary","Results");
	totalNeuriteArea = getResult("Total Area", 0); //get totalNeuriteArea
	fragmentedNeuriteArea = getResult("Total Area", 1); //get fragmentedNeuriteArea
	close("Results");
	
	
	if (saveFragmentationMask) {
		//Save neurite fragmentation masks
		selectWindow("fragment_mask_image");
		run("Cyan");
		run("Invert LUT");
		run("Merge Channels...", "c4=&neurite_mask c6=&fragment_mask_image");
		saveAs("tiff", resultsDir + imageWindow + "_neurite_fragment_mask");
	}
	
	if (saveNucleiMask) {
		//Save nuclei count masks
		selectImage("Drawing of nuclei_mask");
		saveAs("tiff", resultsDir + imageWindow +  "_nuclei_mask");
	}
	
	if (saveNeuriteMask) {
		//Save neurite masks;
		run("Merge Channels...", "c2=[neurite_mask-skeleton] c3=[nuclei_mask] c4=[neurite attachment points mask] create keep");
		selectImage("Composite");
		saveAs("tiff", resultsDir + imageWindow + "_neurite_mask");	
	}
	

	resultString = imageWindow + "," + NeuriteAttachmentPoints+ "," +TotalNeuriteLength + "," + MaxBranchLength+ "," +MeanBranchLength + "," + TotalNeuriteEndPoints+ "," +TotalBranches + "," + TotalTrees+ "," +NucleiCount + "," +NucleiTotalArea + "," + NucleiAverageArea+ "," +TotalNeuronArea+ "," + fragmentedNeuriteArea+ "," +totalNeuriteArea;
	File.append(resultString, resultsFile);
	run("Clear Results");	
	run("Close All");
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


function getImageCalibration(imageWindow) {
	selectWindow(imageWindow);
	getPixelSize(unit, pixelWidth, pixelHeight);
	pixel_size = pixelWidth;
	pixel_unit = unit;
}


function calibrateImage(imageWindow) {
	selectWindow(imageWindow);
	run("Properties...", " unit="+pixel_unit+" pixel_width="+pixel_size+" pixel_height="+pixel_size);
}

function analyzeNeuriteProperties(workingDir) { 
// this function executes the analysis by calling other sub-functions

	experimentId = split(workingDir, "/"); 
	experimentId = experimentId[(lengthOf(experimentId)-1)]; //Extract the experimentId from the directory path
	sourceImagesDir = workingDir + "TIFFs/"; //define source images directory
	fileList = getFileList(sourceImagesDir);
	
	resultsDir = workingDir + "Neurite Analysis_"+ neuriteMarker + "/"; //define output images directory
	resultsFile = workingDir + experimentId + "_" +neuriteMarker + "_Neurite_Analysis.csv";
	
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
	
	
	//Initialise results file
	if (File.exists(resultsFile)) {
		garbage = File.delete(resultsFile); // delete any existing results file
	} 
	f = File.open(resultsFile); //Create a new results file
	print(f, "Filename, NeuriteAttachmentPoints, TotalNeuriteLength, MaxBranchLength, MeanBranchLength, TotalNeuriteEndPoints, TotalBranches, TotalTrees, NucleiCount, NucleiTotalArea, NucleiAverageArea, TotalNeuronArea,FragmentedNeuriteArea,TotalNeuriteArea");
	File.close(f);


	//Analyze images and save results
	for (i = 0; i < fileList.length; i++){
		if (!(endsWith(fileList[i], ".tif"))) {
			exit("Error: Please provide TIFF images for processing"); // return an error and exit the macro if the image is not a TIFF
		}
		filePath = sourceImagesDir + fileList[i];
		open(filePath);
		imageWindow = getTitle();
		if (!(is("composite"))) {
			exit("Error: Please provide composite images for processing"); // return an error and exit the macro if the image is not multichannel
		}
		//Extract the neurite and soma channels from image
		run("Split Channels");
		somaImage = somaChannel + imageWindow;
		neuriteImage = neuriteChannel + imageWindow;
		prepareSomaImage(somaImage);
		prepareNeuriteImage(neuriteImage);
		measureNeuriteProperties(imageWindow, resultsFile, resultsDir);
		showProgress(i, fileList.length);	
	}	
	
	run("Close All");
	print("Neurite Analysis Complete");

}





/*
 	***********************

	Execute code to analyze images and save results

	***********************
*/

analyzeNeuriteProperties(workingDir);













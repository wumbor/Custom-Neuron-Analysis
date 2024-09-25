//This macro measures fluorescence intensity in one channel based on ROIs defined in the other channels of an image hyperstack (i.e. time-lapse image)
//Written by Victor Kumbol, August 2022
//Modified by Victor Kumbol, September 2024



/*
 	***********************

	Set IJ parameters

	***********************
*/
setBatchMode(true);
setOption("ExpandableArrays", true);
run("Clear Results");
run("Set Measurements...", "area mean integrated stack display redirect=None decimal=3");




/*
 	***********************

	Set experiment variables

	***********************
*/
var sizeExclusion = false;
var allChannelColors = newArray("blue", "green", "red", "fred");
var channelOfInterest = "green";
var channelsInImage = 3;
var channelsForCustomROIs =  newArray("red");
allChannelColors = Array.trim(allChannelColors, channelsInImage);





/*
 	***********************

		Functions

	***********************
*/


function AnalyzeKinetics(filePath) {
	//This function performs measurements on all frames in a hyperstack and saves the results to CSV files by calling on other functions
	if (endsWith(filePath, ".nd2")) { //work with both .nd2 files and .tiffs
		run("Bio-Formats Importer", "open=[filePath] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	} else if (endsWith(filePath, ".tif")) {
		open(filePath);
	}

	//Make sure image is a hyperstack
	hyperStack = getTitle();
	if (Stack.isHyperstack) {
		getDimensions(width, height, channels, slices, frames);
	} else {
		exit("The image must be a hyperstack");
	}
	
	//Create an output subfolder
	filePathDir = File.getParent(filePath);
	outputDir = filePathDir + File.separator + File.nameWithoutExtension + "_Results";
	if (!File.exists(outputDir)) {
		File.makeDirectory(outputDir);
		print("");
		print("Output folder created: " + outputDir);
		print("");
	}

	//iterate over all the frames in the hyperstack and perform fluorescence measurements on each
	for (i = 1; i <= frames; i++) {
		selectWindow(hyperStack);
		Stack.setFrame(i);
		run("Reduce Dimensionality...", "channels slices keep");
		currentFrame = getTitle();	
		currentFrameIndex = i;
		MeasureInFrame(hyperStack, currentFrame, currentFrameIndex, outputDir);
		showProgress(i, frames);	
	}

print("Analysis Completed");
print("Total Frames Analysed: " + frames);
}





function MeasureInFrame(hyperStack, frame, frameIndex, outputDir) {
	//This function performs measurements on a single frame and saves the results to a CSV files by calling on other functions

	allChannelNames = newArray("C1-"+frame, "C2-"+frame, "C3-"+frame, "C4-"+frame);
	allChannelNames = Array.trim(allChannelNames, channelsInImage);
	channelOfInterestIndex = getIndexinArray(allChannelColors, channelOfInterest);

	//split the frame into its composite channels
	selectWindow(frame);
	run("Split Channels");
	CustomROIList = CreateCustomROIs(allChannelColors, channelsForCustomROIs); //create ROIs
		
	//Measure fluorescence in all channels using ROI masks from desired channel
	for (i = 0; i < allChannelNames.length; i++) {
		MeasureInROIs(allChannelNames[i], CustomROIList);
		ReformatResults(CustomROIList); //reformat results table for saving
		resultFilePath = outputDir + File.separator + hyperStack +"_miRKinetics_"+ allChannelColors[i]+"-"+frameIndex+".csv";
		selectWindow("Results");
		saveAs("Results", resultFilePath);
	}
	roiManager("Deselect");	
	roiManager("Delete");
	
	selectWindow(hyperStack);
	close("\\Others");
    run("Clear Results");
}








function MeasureInROIs(ChannelName, ROIstoMeasure) { 
	selectWindow(ChannelName); //select the channel to perform measurements on 
	//First, measure fluorescence using individual ROI masks from channels
	for (i = 0; i < lengthOf(ROIstoMeasure); i++) {
		roiManager("select", i);
		//test if ROI name contains "NO_SELECTION" and write NA if true
		name = RoiManager.getName(i);
		if (name.contains("NO_SELECTION_CREATED")) {
			row = nResults;
			setResult("Label", row, "NA");
		    setResult("Area", row, "NA");
		    setResult("Mean", row, "NA");
		    setResult("IntDen", row, "NA");
		    setResult("RawIntDen", row, "NA");
		    setResult("Slice", row, "NA");
		    updateResults();
		} else {
			roiManager("measure");
		}
	}
}


function ReformatResults(roiList) { 
	//Rearrange the results table to display relevant data in desired format
	n = nResults;
  	label = newArray(n);
  	area = newArray(n);
	mean = newArray(n);
	rawintden = newArray(n);
  	slice = newArray(n);

    for (i=0; i<n; i++) {
      label[i] = getResultLabel(i);
      area[i] = getResult("Area", i);
      mean[i] = getResult("Mean", i);
      rawintden[i] = getResult("RawIntDen", i);
      slice[i] = getResult("Slice", i);
    }
   	run("Clear Results");
   	run("Select None");

   for (i=0; i<n; i++) {
      setResult("Channel", i, roiList[i]);
      setResult("Label", i, label[i]);
      setResult("Area", i, area[i]);
      setResult("Mean", i, mean[i]);
      setResult("RawIntDen", i, rawintden[i]);
      setResult("Slice", i, slice[i]);
    }
     updateResults();
}


function CreateROIs(ChannelName, ChannelColor) { 
	ROIs = newArray;
	selectWindow(ChannelName); //select the channel to create masked image
	run("Convert to Mask", "method=Triangle background=Dark calculate");
	if ((sizeExclusion)) {
		run("Analyze Particles...", "size=0-2 show=Masks include in_situ stack"); //Use this to toggle size exclusion
	}
	
	
	 //create ROI of Endosomal Region
	selectWindow(ChannelName); //select the channel to create masked image
	for (i = 1; i <= nSlices; i++) { //loop over the slices, create a selection and ROI from each slice
	    setSlice(i);
	    run("Create Selection");
		roiManager("Add");
				
		type = selectionType(); //test if a valid ROI was created 
		if (type==-1) {
			currentROIindex = (roiManager("count")) - 1;
			roiManager("select", currentROIindex);
			roiManager("rename", "NO_SELECTION_CREATED"); //rename the ROI if invalid
			print(ChannelColor + "-"+ i);
			print("No Selection Created: pHrodo");
		} 			
		ROIs = Array.concat(ROIs, ChannelColor +"-"+ i);
	}
	
	//create ROI of Non-endosomal Region
	selectWindow(ChannelName); //select the channel to create masked image
	for (i = 1; i <= nSlices; i++) { //loop over the slices, create a selection and ROI from each slice
	    setSlice(i);
	    run("Create Selection");
	    run("Make Inverse");
		roiManager("Add");
		
		type = selectionType(); //test if a valid ROI was created 
		if (type==-1) {
			currentROIindex = (roiManager("count")) - 1;
			roiManager("select", currentROIindex);
			roiManager("rename", "NO_SELECTION_CREATED"); //rename the ROI if invalid
			print(ChannelColor + "-"+ i);
			print("No Selection Created: pHrodo Inverse");
		} 			
		ROIs = Array.concat(ROIs, ChannelColor +" inverse -"+ i);
	}
	
	
	selectWindow(ChannelName);
	close();
	return ROIs;
}





function getIndexinArray(array, string) { 
// this function returns the index of a string in an array
	for (i = 0; i < lengthOf(array); i++) {
		if (array[i]==string) {
			index = i;
		}		
	}
	return index;
}


function detectStringinArray(array, string) { 
// this function detects a match of a string in an array
	detected = false;
	for (i = 0; i < lengthOf(array); i++) {
		if (startsWith(string, array[i])) {
			detected = true;
		}		
	}
	return detected;
}


function deleteFolder(folderpath) { 
	//first delete all files in the folder
	deleteFileList = getFileList(folderpath);
		for (i = 0; i < deleteFileList.length; i++){
			deletefilePath = folderpath +deleteFileList[i];
			garbage = File.delete(deletefilePath);
			}	

	//delete the folder itself
	deletecompleted = File.delete(folderpath);
	if (deletecompleted) {
		return true;
	} else {
		return false;
	}
}






 function CreateCustomROIs(allChannelColors, channelsForCustomROIs) {
	//This function creates custom ROIs for measurement, in this case the intersection of the one channel with the complement of another. This practically means the regions where phRodoRed does not intersect with DAPI

	CustomROIList = newArray();
	//Duplicate the channels to be used for creating ROIs
	for (i = 0; i < channelsForCustomROIs.length; i++) { 
	CustomROIChannelIndex = getIndexinArray(allChannelColors, channelsForCustomROIs[i]);
	selectWindow(allChannelNames[CustomROIChannelIndex]);
	run("Duplicate...", "title=StackforROI duplicate"); //Duplicate the channels to be used for creating ROIs
	ROIStackTitle = "StackforROI"; 
	ROIs = CreateROIs(ROIStackTitle, channelsForCustomROIs[i]); //Create ROIs from each stack
	CustomROIList = Array.concat(CustomROIList, ROIs);
	}
	

	return CustomROIList;
}





/*
 	***********************

	Execute code to analyze images and save results

	***********************
*/


//Request file from user and run analysis on all frames in the hyperstack
filePath = File.openDialog("Choose image file"); 
AnalyzeKinetics(filePath);




















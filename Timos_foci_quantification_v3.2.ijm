//choose an input directory of files to be processed
input_dir = getDirectory("Choose input directory");
//choose a directory where ROIs of cells are saved
ROI_dir = getDirectory("Choose a directory of ROIs with the same name as the fluorescent images");
//choose an output directory where the ROI files should be saved
output_dir = getDirectory("Choose output directory");

// percentages to clip
bottomClip = 10;    // clip bottom 10%
topClip = 10;       // clip top 10%
			
// Set the threshold here
Threshold =  1.7;
// Set measurements you want of foci
run("Set Measurements...", "area mean min integrated redirect=None decimal=3");
// set the minimum cell size here (in pixel)
min_cell_size = 50; 

// Empty the ROI manager
wait(500);
nROIs = roiManager("count");
if (nROIs > 1) {
	roiManager("deselect");
	roiManager("Delete");
}
// Close all previously open images
list = getList("image.titles");
if (list.length > 0) {
	for (i = 0; i < list.length; i++) {
		close(list[i]);
	}
}

run("Clear Results");


// Disable antialias correction + padding for binary masks
setOption("BlackBackground", true);
run("Options...", "iterations=0 count=1 black pad do=Nothing");

fileList = getFileList(input_dir);
nFiles = fileList.length;
fileList = getFileList(input_dir);
Array.sort(fileList); // Sort the fileList
nFiles = fileList.length;
ROIList = getFileList(ROI_dir);
Array.sort(ROIList); // Sort the ROIList
nROIFiles = ROIList.length;
print(nFiles + " Images and " + nROIFiles + " sets of ROIs detected.");

// Print the list of files
for (i = 0; i < nFiles; i++) {
    print("Detected file: " + fileList[i]);
    print("Detected file: " + ROIList [i]);
}



for (i = 0; i < nFiles; i ++) { 
	
        file1 = fileList[i];
        fileName = substring(file1, 0, lastIndexOf(file1, ".") );
        file1 = input_dir + fileList[i];
		
        roi = ROIList[i]; 
        roiName = substring(roi, 0, lastIndexOf(roi, ".") );
        roi = ROI_dir + ROIList[i]; 
        
		print(fileName + "   " + roiName);
		
// Check if the file name matches the roi name        
		if (fileName  == roiName) {
        print("Processing " + file1);
        open(file1);
        Image.removeScale();
        fn1=getTitle ();
        rename("img1");
        run("16-bit");
        fileName = substring(fn1, 0, lastIndexOf(fn1, ".") );

		wait(300);
        
// create new image thats black
		selectImage("img1");
        run("Duplicate...", "ignore");
        run("Select All");
		setForegroundColor(0, 0, 0);
		run("Fill", "slice");
		rename("Foci_ROI");
		run("8-bit");
		selectImage("img1");
		
// Get the number of ROIs in the ROI manager
        roiManager("Open", roi);
		nROIs = roiManager("count");

		// Loop through each ROI
		for (j = 0; j < nROIs; j++) {
    		print("Processing ROI index:", j);  // Add this line for debugging
    		
    		 // Remember where the next block of rows will begin
    		startRow = nResults();
			
    		// Select the current ROI
    		selectImage("img1");
    		roiManager("Select", j);
			// duplicate the ROI	
			run("Duplicate...", "title=Cell");
			selectImage("Cell");
			roiManager("Select", j);
			setBackgroundColor(0,0,0);
			run("Clear Outside");
			run("16-bit");
						
			// Get ROI intensity statistics (area, mean, min, max, std, histogram)
			getStatistics(area, mean, min, max, std, histogram); 
			
			// Get total intensity
			IntDen = mean * area;
			
			if (area >= min_cell_size) {
			
			
			total = 0; for (hh = 0; hh < histogram.length; hh++) total += histogram[hh]; // total pixel count (assumed equal to sum(histogram))
			nBins = histogram.length;
			

			// compute bin width and center
			binWidth = (max - min + 1.0) / nBins;
			halfWidth = binWidth / 2.0;
			
			// compute cutoff boundaries in pixel counts
			cutoffLow  = total * (bottomClip / 100.0);                 // skip first bottom %
			cutoffHigh = total - total * (topClip / 100.0);            // upper boundary
			
			accCount = 0.0;
			sum = 0.0;
			includedPixels = 0.0;
			
			for (h = 0; h < nBins; h++) {
			
			     binCount = histogram[h] + 0.0;     // force float
			     binStartCount = accCount;
			     binEndCount   = accCount + binCount;
			
			     // --- compute lowOverlap = min(binEndCount, cutoffHigh)
			     if (binEndCount < cutoffHigh)
			         lowOverlap = binEndCount;
			     else
			         lowOverlap = cutoffHigh;
			
			     // --- compute highOverlap = max(binStartCount, cutoffLow)
			     if (binStartCount > cutoffLow)
			         highOverlap = binStartCount;
			     else
			         highOverlap = cutoffLow;
			
			     // overlap of this bin with keep-range
			     include = lowOverlap - highOverlap;
			
			     // clamp to valid range
			     if (include < 0) include = 0;
			     if (include > binCount) include = binCount;
			
			     // accumulate only if any part of this bin remains
			     if (include > 0) {
			
			         // bin value = center
			         binStartValue = min + h * binWidth;
			         binCenter = binStartValue + halfWidth;
			
			         sum += include * binCenter;
			         includedPixels += include;
			     }
			
			     accCount = binEndCount;
			
			     // early exit when crossing cutoffHigh
			     if (accCount >= cutoffHigh)
			         break;
			}
			
			// guard divide-by-zero
			if (includedPixels > 0)
			     clippedMean = sum / includedPixels;
			else
			     clippedMean = 0;
		
			// Set the threshold and apply to the cell			
			thresh = Threshold * clippedMean;
			
			print("Cell " + j + " has clipped mean " + clippedMean + " threshold of " + thresh);
			selectImage("Cell");
			run("Duplicate...", "title=Cell-1");
			setThreshold(thresh, 65535);
			setOption("BlackBackground", true);
			run("Make Binary");
			run("Create Selection");
			
			if (selectionType != -1) {
				roiManager("Add");
				roiManager("select", nROIs);
				// Only split if the selection is composite
			    if (selectionType == 9) {
			        roiManager("Split");
			        roiManager("select", nROIs);
					roiManager("Delete");
					
					new_nROIs = roiManager("count");
									
					for (h = nROIs; h < new_nROIs; h++) {
						selectImage("Cell");
						roiManager("Select", h);
						run("Measure");
					}
					//waitForUser("wait");
					for (h = new_nROIs -1 ; h >= nROIs; h--) {
						roiManager("Select", h);
						roiManager("Delete");
			    	} 

			    		selectImage("Cell-1");
						close();
			    } else {
			    	selectImage("Cell");
			    	roiManager("select", nROIs);
					run("Measure");
					roiManager("select", nROIs);
					roiManager("delete");
					selectImage("Cell-1");
					close();
					//waitForUser("wait");
				}
			} else {
				r = nResults;
				setResult("Area", r, 0);
				setResult("Mean", r, 0);
				setResult("Min", r, 0);
				setResult("Max", r, 0);
				setResult("IntDen", r, 0);
				setResult("RawIntDen", r, 0);
				//waitForUser("wait");
			}
			
			//waitForUser("wait");
			
    		// Number of rows after this ROI’s processing
    		endRow = nResults();
				
  			 // Add extra metadata columns to *all rows added for this ROI*
		    for (r = startRow; r < endRow; r++) {
		        setResult("Cell_number",    r, "Cell_" + j);
		        setResult("Cell_Area",      r, total);   // 
		        setResult("Cell_Intensity", r, IntDen);  // 
		        setResult("Image", r, fileName);
		        setResult("Clipped_Mean", r, clippedMean);
		        setResult("Threshold", r, thresh);
		    }
				
			// Refresh table display
			updateResults();
			//waitForUser("wait");
			
			// Duplicate foci ROIs to new image
			selectImage("img1");
			run("Select None");
			run("Duplicate...", "ignore");
			setThreshold(thresh, 65535);
			setOption("BlackBackground", true);
			run("Make Binary");
			roiManager("Select", j);
			run("Copy");
			selectImage("Foci_ROI");
			roiManager("select", j);
			setForegroundColor(255, 255, 255);
			run("Paste");
			wait(1000);
			//waitForUser("Foci pasted?");
			selectImage("img1-1");
			wait(500);
			close();

			}
			selectImage("Cell");
			close();
		}


			
			// new image for overlay
			selectImage("img1");
		    run("Duplicate...", "ignore");
		    run("Select All");
			setForegroundColor(0, 0, 0);
			run("Fill", "slice");
			rename("Foci");
			
			selectImage("Foci_ROI");
			run("Create Selection");
			setForegroundColor(80, 80, 80);
			selectImage("Foci");
			run("Restore Selection");
			run("Draw", "slice");
			run("Fill", "slice");			
			run("16-bit");
			
			run("Merge Channels...", "c4=img1 c6=Foci create keep");
			run("Channels Tool...");
			// waitForUser("Are you happy with the segmentation?");
			save(output_dir + fileName + "_foci_overlay.tif");
			selectImage("Foci_ROI");
			save(output_dir + fileName + "_foci.tif");
			close("*");
			roiManager("deselect");
			roiManager("delete");
		
	}
}		

saveAs("Results", output_dir + "Results.csv");






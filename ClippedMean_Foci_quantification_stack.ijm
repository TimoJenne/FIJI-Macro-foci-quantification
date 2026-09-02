//choose an input directory of files to be processed
input_dir = getDirectory("Choose input directory");
//choose a directory where ROIs of cells are save
ROI_dir = getDirectory("Choose a directory of ROIs with the same name as the fluorescent images");
//choose an output directory where the ROI files should be saved
output_dir = getDirectory("Choose output directory");

// percentages to clip
bottomClip = 35;    // clip bottom 10%
topClip = 10;       // clip top 10%			
// Set the threshold here
Threshold =  2.4;

//setBatchMode(true); 


// Set measurements you want of foci
run("Set Measurements...", "area mean min integrated redirect=None decimal=3");
// set the minimum cell size here (in pixel)
min_cell_size = 1000; 

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



if (nFiles == nROIFiles ) {
	for (i = 0; i < nFiles; i++) {
        file1 = input_dir + fileList[i];

        roi = ROI_dir + ROIList[i ]; 
        
        if (File.exists(file1) && File.exists(roi)) {

        fileName = substring(file1, 0, lastIndexOf(file1, ".") );     

        roiName = substring(roi, 0, lastIndexOf(roi, ".") );
        
		print(fileName + "   " + roiName);
		   

        print("Processing " + file1);
        open(file1);
        fn1=getTitle ();
        rename("img1");
        
        
        run("Split Channels");
		selectImage("C1-img1");
		rename("img1");
		Image.removeScale();
        run("16-bit");
		
		selectImage("C2-img1");
		rename("img2");
		Image.removeScale();
        run("16-bit");
        
        fileName = substring(fn1, 0, lastIndexOf(fn1, ".") );
                
        //waitForUser("wait");
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
    								
			// Get ROI intensity statistics (area, mean, min, max, std, histogram)
			getStatistics(area, mean, min, max, std, histogram); 

			if (area >= min_cell_size) {
				
			// duplicate the ROI	
			run("Duplicate...", "title=Cell");
			selectImage("Cell");
			roiManager("Select", j);
			setBackgroundColor(0,0,0);
			run("Clear Outside");
			run("16-bit");
			
			// duplicate the cell in channel 2 as well
			selectImage("img2");
    		roiManager("Select", j);
    		run("Duplicate...", "title=Cell_Ch2");
			selectImage("Cell_Ch2");
			roiManager("Select", j);
			setBackgroundColor(0,0,0);
			run("Clear Outside");
			run("16-bit");
			Ch2Int = getValue("IntDen");
    		
			
			// Get total intensity
			IntDen = mean * area;
			
			
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
			if (includedPixels > 0) {
			     clippedMean = sum / includedPixels;
			} else {
			     clippedMean = 0;
			}
			// Set the threshold and apply to the cell			
			thresh = Threshold * clippedMean;
			
			print("Cell " + j + " has clipped mean " + clippedMean + " threshold of " + thresh);
			run("Select None");
			wait(100);
			selectImage("Cell");
			run("Duplicate...", "title=Cell-1");
			
			cellDupTitle = getTitle();
			if (!isOpen("Cell-1") || cellDupTitle == "Cell") {
			    print("WARNING: Duplicate failed for ROI " + j + " in " + fileName + " — skipping this cell.");
			    r = nResults;
			    setResult("Area", r, -1);
			    setResult("Cell_number", r, "Cell_" + j);
			    setResult("Image", r, fileName);
			    selectImage("Cell");
			    close();
			    continue;
			}

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
						//
						endRow = nResults();
					    for (r = startRow; r < endRow; r++) {
					        setResult("Cell_number",    r, "Cell_" + j);
					        setResult("Cell_Area",      r, total);   // 
					        setResult("Cell_Intensity", r, IntDen);  // 
					        setResult("Image", r, fileName);
					        setResult("Clipped_Mean", r, clippedMean);
					        setResult("Threshold", r, thresh);
					        setResult("Channel", r, 1);
					    }
					updateResults();
					wait(200);
					//waitForUser("wait");
					
					//measure intensity in channel 2
					selectImage("Cell_Ch2"); 	
					startRow = nResults;
					
					for (h = nROIs; h < new_nROIs; h++) {
						selectImage("Cell_Ch2");
						roiManager("Select", h);
						run("Measure");
						}
					
					endRow = nResults();
					    for (r = startRow; r < endRow; r++) {
					        setResult("Cell_number",    r, "Cell_" + j);
					        setResult("Cell_Area",      r, total);   // 
					        setResult("Cell_Intensity", r, Ch2Int);  // 
					        setResult("Image", r, fileName);
					        setResult("Clipped_Mean", r, clippedMean);
					        setResult("Threshold", r, thresh);
					        setResult("Channel", r, 2);
					    }
					    
					updateResults();
					wait(200);
					//waitForUser("wait");
					for (h = new_nROIs -1 ; h >= nROIs; h--) {
						roiManager("Select", h);
						roiManager("Delete");
			    	} 

			    		selectImage("Cell-1");
						close();
						selectImage("Cell_Ch2");
						close();
			    } else {
			    	selectImage("Cell");
			    	roiManager("select", nROIs);
					run("Measure");
					 
					endRow= nResults;
					//
					for (r = startRow; r < endRow; r++) {
				        setResult("Cell_number",    r, "Cell_" + j);
				        setResult("Cell_Area",      r, total);   // 
				        setResult("Cell_Intensity", r, IntDen);  // 
				        setResult("Image", r, fileName);
				        setResult("Clipped_Mean", r, clippedMean);
				        setResult("Threshold", r, thresh);
				        setResult("Channel", r, 1);
				    }
				    updateResults();
				    wait(100);
				    
				    startRow = nResults;
					selectImage("Cell_Ch2");
			    	roiManager("select", nROIs);
					run("Measure");
					 
					endRow= nResults;
					//
					for (r = startRow; r < endRow; r++) {
				        setResult("Cell_number",    r, "Cell_" + j);
				        setResult("Cell_Area",      r, total);   // 
				        setResult("Cell_Intensity", r, Ch2Int);  // 
				        setResult("Image", r, fileName);
				        setResult("Clipped_Mean", r, clippedMean);
				        setResult("Threshold", r, thresh);
				        setResult("Channel", r, 2);
				    }					
				    updateResults();
				    wait(100);					
					
					roiManager("select", nROIs);
					roiManager("delete");
					selectImage("Cell-1");
					close();
					
					selectImage("Cell_Ch2");
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
				setResult("Cell_number",    r, "Cell_" + j);
				setResult("Cell_Area",      r, total);   // 
				setResult("Cell_Intensity", r, IntDen);  // 
				setResult("Image", r, fileName);
				setResult("Clipped_Mean", r, clippedMean);
				setResult("Threshold", r, thresh);
				setResult("Channel", r, 1);
				
				r = nResults;
				setResult("Area", r, 0);
				setResult("Mean", r, 0);
				setResult("Min", r, 0);
				setResult("Max", r, 0);
				setResult("IntDen", r, 0);
				setResult("RawIntDen", r, 0);
				setResult("Cell_number",    r, "Cell_" + j);
				setResult("Cell_Area",      r, total);   // 
				setResult("Cell_Intensity", r, Ch2Int);  // 
				setResult("Image", r, fileName);
				setResult("Clipped_Mean", r, clippedMean);
				setResult("Threshold", r, thresh);
				setResult("Channel", r, 2);
				
				
				selectImage("Cell-1");
				close();
				selectImage("Cell_Ch2");
				close();
				//waitForUser("wait");
			}
			updateResults();
			wait(200);

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
			
						selectImage("Cell");
			close();
			}

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
			
			run("Merge Channels...", "c4=img1 c6=Foci 	create 	keep");
			//run("Channels Tool...");
			//waitForUser("Are you happy with the segmentation?");
			save(output_dir + fileName + "_foci_overlay.tif");
			
			run("Merge Channels...", "c4=img2 c6=Foci create keep");
			save(output_dir + fileName2 + "_foci_overlay.tif");
			
			selectImage("Foci_ROI");
			save(output_dir + fileName + "_foci.tif");
			close("*");
			roiManager("deselect");
			roiManager("delete");
        } else {
			print("Files don't exist");
        }
	}
}	else {
print("Number of files and ROIs don't match");	
}

saveAs("Results", output_dir + "Results.csv");






// Spot finder script
// -------------------------------
// 
// This script takes a 3D+time+channel stack as input. It determines bright 
// spots in the centriole_channel (configured below) and writes the coordinates 
// to a tableassumes.
// 
// Installation: 
//   To use this script in Fiji (http://fiji.sc/Downloads), you need to activate
//   the clij and clij2 update sites in the menu Help > Update... under 
//   Manage Update Sites. Restart Fiji after installation. 
//   Read more about CLIJ online: http://clij.github.io
//
// Usage:
//   * Drag and drop this file on the Fiji main window.
//   * Configure filename and properts below this comment.
//   * Click "run" in Fijis script editor.
//
// Author: Robert Haase, rhaase@mpi-cbg.de.
// April 2020
// License: BSD3
// 
// Copyright 2020 Robert Haase, Max Planck Institute for Molecular Cell Biology and Genetics Dresden
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//Close and reset all windows
close("*");
roiManager("reset");
run("Clear Results");
//
// configure input file. Do not use \ . Use / instead. 
filename =         "D:/data/Analysis Data/Exp_data/Exp_RNAi_c27d9.1/Exp59a_tmr31-c27d9-1-08cell_LLSM_20220126/Exp59a_tmr31-c27d9-1-08cell_LLSM_20220126_Analysis/8cell_ABar/8cell_ABar_ImageJ/8cell_ABar.tif"

// configure where results should be saved
coordinates_csv =  "D:/data/Analysis Data/Exp_data/Exp_RNAi_c27d9.1/Exp59a_tmr31-c27d9-1-08cell_LLSM_20220126/Exp59a_tmr31-c27d9-1-08cell_LLSM_20220126_Analysis/8cell_ABar/8cell_ABar_ImageJ/8cell_ABar_coordinates.csv";
coordinates_rois = "D:/data/Analysis Data/Exp_data/Exp_RNAi_c27d9.1/Exp59a_tmr31-c27d9-1-08cell_LLSM_20220126/Exp59a_tmr31-c27d9-1-08cell_LLSM_20220126_Analysis/8cell_ABar/8cell_ABar_ImageJ/8cell_ABar_coordinates_roi.zip";

// configure a folder to output a video; must end with /
coordinates_video = "D:/data/Analysis Data/Exp_data/Exp_RNAi_c27d9.1/Exp59a_tmr31-c27d9-1-08cell_LLSM_20220126/Exp59a_tmr31-c27d9-1-08cell_LLSM_20220126_Analysis/8cell_ABar/8cell_ABar_ImageJ/8cell_ABar_video/"

// channel configuration, first channel is 1
centriole_channel = 2; 
dna_channel = 1;

// enter correct pixel size here
pixel_size_x = 0.1;
pixel_size_y = 0.1;
pixel_size_z = 0.5;

zoom = 20; // one micron is divided into <zoom> pixels

// clean up at the beginning
run("Close All");

// Init GPU
run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();

// init ROI Manager
roiManager("reset");

// initiate coordinates table on disc
if (File.exists(coordinates_csv)) {
	File.delete(coordinates_csv);
}
File.append("X1,Y1,Z1,X2,Y2,Z2, frame, minimum distance index, second distance minimum index", coordinates_csv);

// initiate video folder
File.makeDirectory(coordinates_video); 

// load data
open(filename);
run("Select None");

// read image properties
input = getTitle();
getDimensions(width, height, channels, slices, frames);
center_x = width / 2;
center_y = height / 2;
center_z = slices / 2;
projections = newArray("x", "y", "z");



for (frame = 0; frame < frames; frame++) {
	
	//print the frame number on the log window
	print(frame, "frame");
	
	// select frame and channel; push stack to GPU
	selectWindow(input);
	Stack.setFrame(frame + 1);	
	Stack.setChannel(centriole_channel);
	Ext.CLIJ2_pushCurrentZStack(input);
	Ext.CLIJ2_pull(input);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	close();
	print("Statistics of the current image", mean, min, max, stdDev); //stdD
	//Ext.CLIJ2_getMeanOfAllPixels(input, mean);

	//Detect the intensity of the spot by multiplying the stdev by a factor of 2 or more. The higher the better.
	threshold = mean + 1.3*stdDev;
	//threshold = mean + 2*stdDev;
	//threshold = mean + 2.5*stdDev;

	// gaussian blur
	blurred = "blurred";
//	Ext.CLIJ2_gaussianBlur3D(input, blurred, 1.1, 1.1, 0);
	Ext.CLIJ2_gaussianBlur3D(input, blurred, 1.3, 1.3, 0);
	//Ext.CLIJ2_median3DSphere(input, blurred, 1.1, 1.1, 0);
	//Ext.CLIJ2_median3DSphere(input, blurred, 2, 2, 0);
	//Ext.CLIJ2_median3DSphere(input, blurred, 4, 4, 0);

	// thresholding
//	thresholded = "thresholded";
	//Ext.CLIJ2_thresholdMaxEntropy(blurred, thresholded);
	//Ext.CLIJ2_greaterConstant(blurred, thresholded, threshold);
	Ext.CLIJ2_automaticThreshold(blurred, thresholded, "Phansalkar");
	//Ext.CLIJ2_thresholdMean(blurred, thresholded);

	// maxima detection
	maxima = "maxima";
	radius = 5.0; // !check the radius adjustment to make detection more robust
	Ext.CLIJ2_detectMaximaBox(blurred, maxima, radius);
	
//	Ext.CLIJ2_detectMaxima3DBox(blurred, maxima, radius, radius, radius);

	// eliminate maxima out of the thresholded region (noise)
	masked_maxima = "masked_maxima";
	Ext.CLIJ2_mask(maxima, thresholded, masked_maxima);


	// get list of detected points into the results table
	pointlist = "pointlist";
	Ext.CLIJ2_spotsToPointList(masked_maxima, pointlist);
	pointlist_xyz = "pointlist_xyz";
	Ext.CLIJ2_transposeXY(pointlist, pointlist_xyz);
	resultsBefore = nResults(); 
	Ext.CLIJ2_image2DToResultsTable(pointlist_xyz);
	resultsAfter = nResults();
	Ext.CLIJ2_release(pointlist);
	Ext.CLIJ2_release(pointlist_xyz);
	print("Result_Before_After", resultsAfter, resultsBefore);

	// determine the two points which are the closest to the image stack center
	minimum_distance = 0;
	minimum_distance_index = -1;
	second_minimum_distance = 0;
	second_minimum_distance_index = -1;
	
	if (frame == 0) {

		for (i = resultsBefore; i < resultsAfter; i++) {
			x = getResult("X0", i);
			y = getResult("X1", i);
			z = getResult("X2", i);
	
			// determine distance from center to eliminate distant wrong detections
		    //copy this code from 157 to 167  
		    distance_from_center = sqrt(
		    	pow(pixel_size_x * (center_x - x), 2) + 
		    	pow(pixel_size_y * (center_y - y), 2) + 
		    	pow(pixel_size_z * (center_z - z), 2)
		    );
	
		    //save the distance_from_center on the result table
		    setResult("Distance_from_center", i, distance_from_center);
			setResult("Frame", i, frame);
			
			//deteremine the minimum distance between the two detected points on the first frame
			if (distance_from_center < minimum_distance || minimum_distance_index == -1) {
				second_minimum_distance = minimum_distance;
				second_minimum_distance_index = minimum_distance_index;
				minimum_distance = distance_from_center;
				minimum_distance_index = i;
			} else if (distance_from_center < second_minimum_distance || second_minimum_distance_index == -1) {
				second_minimum_distance = distance_from_center;
				second_minimum_distance_index = i;
			}
		}
		
		//get the x,y,z coordinates from the result table
		x1 = getResult("X0", minimum_distance_index);
		y1 = getResult("X1", minimum_distance_index);
		z1 = getResult("X2", minimum_distance_index);
		x2 = getResult("X0", second_minimum_distance_index);
		y2 = getResult("X1", second_minimum_distance_index);
		z2 = getResult("X2", second_minimum_distance_index);
		
		//show the coordinates of the 2 detected points at the x,y,z 			
		print("point1: " + x1 + "/" + y1 + "/ " + z1);
		print("point2: " + x2 + "/" + y2 + "/ " + z2);

		//set the position of the points on the 4D stack and add on the ROI Manager
		//first point
		makePoint(x1, y1);
		Roi.setPosition(2, z1 + 1, 1);
		roiManager("add");
	
		//second point
		makePoint(x2, y2);
		Roi.setPosition(2, z2 + 1, 1);
		roiManager("add");
	
		//ask for validation for correct point detection
		//if the detected spot is incorrect, move the incorrect spot to the correct position and update on the roi manager window
		waitForUser("Validation", "Is the detected spots correct?\n\nPlease click the point in the roiManager.\n\nIf the detected spot is incorrect, move the incorrect spot to the correct position and update on the roi manager window.\n\nClick OK afterward");
				
		//select the coordinate of the first point from the roi manager after validation
		roiManager("select", 0);
		Roi.getCoordinates(xpoints, ypoints);
		x1 = xpoints[0];
		y1 = ypoints[0];
		Stack.getPosition(_, z1, _);
		//subtract to get the correct position of the z on the roiManger
		z1 = z1 - 1;
				
		//select the coordinate of the second point from the roi manager after validation
		roiManager("select", 1);
		Roi.getCoordinates(xpoints, ypoints);
		x2 = xpoints[0];
		y2 = ypoints[0];
		Stack.getPosition(_, z2, _);
		//subtract to get the correct position of the z on the roiManger
		z2 = z2 - 1;

		//save the x,y,z coordinates on the result table
		setResult("X0", minimum_distance_index, x1);
		setResult("X1", minimum_distance_index, y1);
		setResult("X2", minimum_distance_index, z1);
		setResult("X0", second_minimum_distance_index, x2);
		setResult("X1", second_minimum_distance_index, y2);
		setResult("X2", second_minimum_distance_index, z2);
		
		//show the coordinates of the 2 detected points at the x,y,z
		print("point1: " + x1 + "/" + y1 + "/ " + z1);
		print("point2: " + x2 + "/" + y2 + "/ " + z2);
						
	}	else {
		//set the condition for the detecting the points in the other frames after validating the detection on the first frame
		for (i = resultsBefore; i < resultsAfter; i++) {

			x = getResult("X0", i);
			y = getResult("X1", i);
			z = getResult("X2", i);
	
			// determine distance from first_point to eliminate distant wrong detections
		    distance_from_first_point = sqrt(
		    	pow(pixel_size_x * (former_x1 - x), 2) + 
		    	pow(pixel_size_y * (former_y1 - y), 2) + 
		    	pow(pixel_size_z * (former_z1 - z), 2)
		    );
		    
			distance_from_second_point = sqrt(
		    	pow(pixel_size_x * (former_x2 - x), 2) + 
		    	pow(pixel_size_y * (former_y2 - y), 2) + 
		    	pow(pixel_size_z * (former_z2 - z), 2)
		    );
		    //save the Distance_from_first_point on the result table
		    setResult("Distance_from_first_point", i, distance_from_first_point);
			setResult("Frame", i, frame);

			//save the Distance_from_second_point on the result table
			setResult("Distance_from_second_point", i, distance_from_second_point);
			
			//deteremine the minimum distance between the two detected points on the other frames
			if (distance_from_first_point < minimum_distance || minimum_distance_index == -1) {
				minimum_distance = distance_from_first_point;
				minimum_distance_index = i;
			} 
			if (distance_from_second_point < second_minimum_distance || second_minimum_distance_index == -1) {
				second_minimum_distance = distance_from_second_point;
				second_minimum_distance_index = i;
			}
		}
	}
		
	if (minimum_distance_index == -1 || second_minimum_distance_index == -1) {
	
		//In this case, if point detection didn't work
		print("Point detection not found in frame" + frame);
		File.append(0 + "," + 0 + "," + 0  + "," + 0 + "," + 0 + "," + 0 + ","+frame+"," + minimum_distance_index + "," + second_minimum_distance_index, coordinates_csv);
	} else {

		setResult("Distance_order", minimum_distance_index, 1);
		setResult("Distance_order", second_minimum_distance_index, 2);
		print("Line indexes", minimum_distance_index, second_minimum_distance_index);
	
		// the two points are:
		x1 = getResult("X0", minimum_distance_index);
		y1 = getResult("X1", minimum_distance_index);
		z1 = getResult("X2", minimum_distance_index);
		x2 = getResult("X0", second_minimum_distance_index);
		y2 = getResult("X1", second_minimum_distance_index);
		z2 = getResult("X2", second_minimum_distance_index);			
		print("point1: " + x1 + "/" + y1 + "/ " + z1);
		print("point2: " + x2 + "/" + y2 + "/ " + z2);
	 	
		//redefine the coordinates of the two points
		former_x1 = x1;
		former_y1 = y1;
		former_z1 = z1;
		former_x2 = x2;
		former_y2 = y2;
		former_z2 = z2;
	
		makePoint(x1, y1);
		Roi.setPosition(2, z1 + 1, frame + 1);
		roiManager("add");
		roiManager("select", roiManager("count") - 1);
		roiManager("rename", "First_Spot" + frame);
	
		// save video
		run("Flatten", "slice");
		number = "000000" + (frame * 2);
		number = substring(number, lengthOf(number) - 6);
		save(coordinates_video + number + ".tif");
		close();
	
		makePoint(x2, y2);
		Roi.setPosition(2, z2 + 1, frame + 1);
		roiManager("add");
		roiManager("select", roiManager("count") - 1);
		roiManager("rename", "Second_Spot" + frame);
		
	
		// save video
		run("Flatten", "slice");
		number = "000000" + (frame * 2 + 1);
		number = substring(number, lengthOf(number) - 6);
		save(coordinates_video + number + ".tif");
		close();
		
		
		File.append(x1 + "," + y1 + "," + z1  + "," + x2 + "," + y2 + "," + z2 + ","+frame+"," + minimum_distance_index + "," + second_minimum_distance_index, coordinates_csv);
		
	} 
	if (minimum_distance_index == second_minimum_distance_index) {
		print("Error: Both point are the same");
		waitForUser("Error: Both point are the same");
		exit
	}
}

//save roi manager
roiManager("save", coordinates_rois);

// clean up and say bye
Ext.CLIJ2_reportMemory();
Ext.CLIJ2_clear();
print("Bye");



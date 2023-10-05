//Kymograph generation script

//clean up and reset all windows
close("*");
roiManager("reset");

//define the path
folder = "D:/data/Analysis Data/Exp_data/Exp_RNAi_c27d9.1/Exp59a_tmr31-c27d9-1-08cell_LLSM_20220126/Exp59a_tmr31-c27d9-1-08cell_LLSM_20220126_Analysis/8cell_ABar/8cell_ABar_ImageJ/";
//image = "8cell_ABar_kymo/8cell_ABar.tif";
image = "8cell_ABar.tif";
table = "8cell_ABar_coordinates.csv";
Kymo = "8cell_ABar_kymo/8cell_ABar_kymo.tif";

//loading the image
open(folder + image);
open(folder + table);

//open(image);
//open(table);

//open(table);
Table.rename(table, "Results");

// one micron is divided into pixels

//Init GPU 
run("CLIJ2 Macro Extensions", "cl_device="); 
Ext.CLIJ2_clear();

//read image properties 
input = getTitle(); 
getDimensions(width, height, channels, slices, frames);
center_x = width / 2; 
center_y = height / 2; 
center_z = slices / 2;
projections = newArray("x", "y", "z");

// channel configuration, first channel is 1
centriole_channel = 2; 
dna_channel = 1;

// enter correct pixel size here
pixel_size_x = 0.1;
pixel_size_y = 0.1;
pixel_size_z = 0.5;
zoom = 10;

// read through the frames for each file
for (frame = 0; frame < frames; frame++) {

	//select frame and channel; push stack to GPU
	selectWindow(input);    
	Stack.setFrame(frame + 1);
	
	//Stack.setChannel(centriole_channel);
	Ext.CLIJ2_pushCurrentZStack(input);

	//print out the coordinates from the table, convert the coordinates to integers and rename the columns
	x1 = parseInt(getResult('X1', frame)); 
	y1 = parseInt(getResult('Y1', frame));
	z1 = parseInt(getResult('Z1', frame));
	x2 = parseInt(getResult('X2', frame));
	y2 = parseInt(getResult('Y2', frame));
	z2 = parseInt(getResult('Z2', frame));
	
	//add to roiManager for both points
	makePoint(x1, y1);
	Roi.setPosition(2, z1 + 1, frame +1); // 1 has to be replaced by a frame number (HINT)
	roiManager("add");
	
	makePoint(x2, y2);
	Roi.setPosition(2, z2 + 1, frame +1); // 1 has to be replaced by a frame number (HINT)
	roiManager("add");

	// determine the center of the detected line. 
	// This will be the center of the isotropic transformed stack
	mid_x = (x1 + x2) / 2;
	mid_y = (y1 + y2) / 2;
	mid_z = (z1 + z2) / 2;
	
	// determine angles to rotate the stack
	deltaX = (x2 - x1) * pixel_size_x;
	deltaY = (y2 - y1) * pixel_size_y;
	deltaZ = (z2 - z1) * pixel_size_z;
	distXY = sqrt(pow(deltaX, 2) + pow(deltaY, 2));
	distXZ = sqrt(pow(deltaX, 2) + pow(deltaZ, 2));
	distYZ = sqrt(pow(deltaY, 2) + pow(deltaZ, 2));

	// length
	length = sqrt(pow(deltaX, 2) + pow(deltaY, 2) + pow(deltaZ, 2));
	
	// rotate in xy plane around z
	angleXY = -Math.atan2(deltaY, deltaX) / PI * 180;
	
	// rotate xz plane around y
	// angleXZ = -Math.atan2(sqrt(pow(deltaX, 2) + pow(deltaY, 2)), pow(deltaZ, 2)) / PI * 180;
	angleXZ = -Math.acos(deltaZ / length) / PI * 180; 

	// keep the x projection constant
	angleYZ = 0;
	
	// determine size of the new isotropic stack    
	newWidth = width * pixel_size_x * zoom;
	newHeight = height * pixel_size_y * zoom;
	
	// elongate stack in Z to not miss centrioles moving apart (cropping a bit wider before and after the 2 points to avoid missing the poles)
	elongation_z = 2.0;
	//elongation_z = 4.0;
	newDepth = slices * pixel_size_z * zoom * elongation_z; 
//	print(newDepth);
	
	scale_factor_x = newWidth / width;
	scale_factor_y = newHeight / height;
	scale_factor_z = newDepth / slices / elongation_z;
	
	
//	print(scale_factor_z);
	
	// formulate the affine transform
	// see https://clij.github.io/clij2-docs/reference_affineTransform3D
	transform = 
		// scale the stack to isotropic voxels
		" scaleX=" + (scale_factor_x) +
		" scaleY=" + (scale_factor_y) +
		" scaleZ=" + (scale_factor_z) + 
		// translate it to the center point between the two detect points
		" translateX=" + (mid_x * scale_factor_x) + 
		" translateY=" + (mid_y * scale_factor_y) + 
		" translateZ=" + (mid_z * scale_factor_z) + 
		// rotate it around the center
		" rotateZ=" + (angleXY) + 
		" rotateY=" + (angleXZ) + 
		" rotateX=" + (angleYZ) + 
		
		// move it back to the center of the new image stack
		" translateX=" + (-center_x * scale_factor_x) + 
		" translateY=" + (-center_y * scale_factor_y) + 
		" translateZ=" + (-center_z * scale_factor_z * elongation_z) + 
		"";
	
	
	// apply transform and maximum projections to all channels
	for (c = 0; c < channels; c++) {
		
		// push specific channel stack to GPU
		selectWindow(input);
		Stack.setChannel(c + 1);
		Ext.CLIJ2_pushCurrentZStack(input);
		
		// actually apply transform
		transformed_input = "transformed_input";
		Ext.CLIJ2_create3D(transformed_input, newWidth, newHeight, newDepth, 32);   
		Ext.CLIJ2_affineTransform3D(input, transformed_input, transform);
		print(transform);
		Ext.CLIJ2_transposeXZ(transformed_input, rotated_input);
		
		// crop out the center part of the stack
		cropped_stack = "cropped_stack";
		
		//Ext.CLIJ2_crop3D(transformed_input, cropped_stack, newWidth / 4, newHeight / 4, 0, newWidth / 2, newHeight / 2, newDepth); 
		//Ext.CLIJ2_crop3D(transformed_input, cropped_stack, newWidth / 6*2, newHeight / 6*2, 0, newWidth / 3, newHeight / 3, newDepth); 
		//Ext.CLIJ2_crop3D(transformed_input, cropped_stack, newWidth / 8*3, newHeight / 8*3, 0, newWidth / 4, newHeight / 4, newDepth);
		Ext.CLIJ2_crop3D(transformed_input, cropped_stack, newWidth / 16*7, newHeight / 16*7, 0, newWidth / 8, newHeight / 8, newDepth);
		//Ext.CLIJ2_crop3D(transformed_input, cropped_stack, newWidth / 32*9, newHeight / 32*9, 0, newWidth / 16, newHeight / 16, newDepth);
		//Ext.CLIJ2_pull(cropped_stack);
		
		// measure intensity per plane
		measurements = "measurements";
		Ext.CLIJ2_sumImageSliceBySlice(cropped_stack, measurements);	

		// save measurments in a kymograph
		measurement_kymograph = "measurement_kymograph" + c;
		Ext.CLIJ2_create3D(measurement_kymograph, newDepth, frames, 1, 32);
		Ext.CLIJ2_paste3D(measurements, measurement_kymograph, 0, frame, 0);
	
	// generate maximum X, Y and Z projections of the isotropic stack
		for (p = 0; p < 3; p++) {
		projection = projections[p];
		max_projection = "max_projection_" + projection;
		max_projection_video = "max_projection_video_" + projection + c;
		
			if (p == 0) {
			    Ext.CLIJ2_maximumXProjection(transformed_input, max_projection);
			    Ext.CLIJ2_create3D(max_projection_video, newDepth, newHeight, frames, 32);
			}else if (p == 1) {
			    Ext.CLIJ2_maximumYProjection(transformed_input, max_projection);
			    Ext.CLIJ2_create3D(max_projection_video, newWidth, newDepth, frames, 32);
			}else if (p == 2) {
			    Ext.CLIJ2_maximumZProjection(transformed_input, max_projection);
			    Ext.CLIJ2_create3D(max_projection_video, newWidth, newHeight, frames, 32);
			}
			
				Ext.CLIJ2_copySlice(max_projection, max_projection_video, frame);
		}
	}
}

	
//show and merge kymographs 
merge_string = "" 
for (c = 0; c < channels; c++) { 
	measurement_kymograph = "measurement_kymograph" + c; 
	Ext.CLIJ2_pull(measurement_kymograph); 
	merge_string = merge_string + "c" + (c+1) + "=" + measurement_kymograph + " ";
	 
} 
merge_string = merge_string + "create"; 
run("Merge Channels...", merge_string); 
rename("kymograph");
Stack.setChannel(1);  
run("Enhance Contrast", "saturated=0.35");
Stack.setChannel(2);
run("Enhance Contrast", "saturated=0.35");
saveAs("tif", folder + Kymo);
//saveAs("tif", Kymo);

// show and merge maximum projections 
for (p = 0; p < 3; p++) { 
	merge_string = ""; 
	projection = projections[p]; 
	for (c = 0; c < channels; c++) { 
		max_projection_video = "max_projection_video_" + projection + c; 
		Ext.CLIJ2_pull(max_projection_video); 
		merge_string = merge_string + "c" + (c+1) + "=" + max_projection_video + " ";
		 
	} 
	merge_string = merge_string + "create"; 
	run("Merge Channels...", merge_string); 
	rename("Max " + projection + " projection");	
	Stack.setChannel(1);
	run("Enhance Contrast", "saturated=0.35"); 
	Stack.setChannel(2);
	run("Enhance Contrast", "saturated=0.35");
}
Ext.CLIJ2_reportMemory(); 
Ext.CLIJ2_clear(); 
print("Bye");



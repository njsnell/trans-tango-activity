// Registered_BNW_to_Heatmap

//input = getDirectory("Input directory");
//output = getDirectory("Output directory");
//bgImage_path = File.openDialog("Choose background image");

base_input = getDirectory("Choose fly folder");
input = base_input + "Heatmaps/";
output = base_input + "Images/";

processFolder(input);

close("*");

function processFolder(input) {
	list = getFileList(input);
	for (k = 0; k < list.length; k++) {
		if(endsWith(list[k], "/"))
			processFolder("" + input + list[k]);
		else if(endsWith(list[k], ".tif"))
			processTIF(input, list[k]);
	}
}

function processTIF(input, file) {
	open(input + file);
	run("physics");
	setMinAndMax(0,300);
	png_filename = replace(file,".tif",".png");
	saveAs("PNG", output + png_filename);
	close();
}

// Choose input and output folders
// Choose background image
//
// For each tif in input folder:
//		Open tif
//		Merge tif with background image
//		Run physics LUT
//		Set max value
//		Save composite as .png in output folder
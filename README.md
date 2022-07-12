# trans-tango-activity

Code used for data analysis in "Complex Representation of Taste Quality by Second-Order Gustatory Neurons in Drosophila" by Snell et al.

Contains MATLAB scripts used for analysis (A1,A2,A3) plus associated MATLAB functions, and ImageJ macro for heatmap visualization

Three main MATLAB scripts:

1) A1_Preprocessing

Used to preprocess imaging data.

Input: Folder with tif stacks of imaging data, each stack containing 1 individual trial from 1 fly

Outputs: mocoPlanesAll.mat file (.mat file containing all motion-corrected imaging data from 1 fly), SC (spatial components of ROIs from 1 fly), TC (temporal components of ROIs from 1 fly)

2) A2_AnalyzeROIs

Used for ROI-based analysis of imaging data.

Input: Cell arrays of SC (spatial components) and TC (temporal components), the outputs from A1 script.

3) A3_Heatmaps

Used to generate heatmaps to visualize stimulus response.

Input: mocoPlanesAll.mat file generated in A1 script.
 

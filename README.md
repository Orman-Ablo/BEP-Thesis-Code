# BEP-Thesis-Code
A collection of Code used to Produce the OTF fitting Routine used in the BEP thesis

In order to use The second fit, Follow the pipeline set for the 'vectorfit' Project:
https://gitlab.tudelft.nl/imphys/ci/vectorfit

Alternatively, the code for the 'Vectorfit' Project can be found in the 'vectorfit-master'.zip folder.


Once you have collected the initial fit of your localizations, follow the Pipeline for the OTF fit:

1. run the file 'combining_locations.m' or 'combined_locations.m';, which isolate the localizations of all ROIs. 
2. run 'PSF_data.m': it produces the 3-D stacks of the PSFs required for the 3-D OTF. you can change the upscaling factor to  produce finer/rougher sampling of the PSF. You need to add in the combined locations first.
3. run 'Plot_PSF_Results.m': this file produces the used 3-D OTF, and plots the shapes of the PSFs and OTFs. You can define here the desired size of the OTF
4. run 'OTF_Fit_All.m': it is a modified copy of the 'localize.m'file used in the initial fit. Set params.mult_run and params.FlagOTF as true to run the OTF fit.
If you want to produce any of the plots used in the report, their pipeline is listed below:

CRLB:
1. run 'CRLB_calc': it combines the CRLB, and the positions
2. run 'plot_CRLB': Produces the figure of the mean CRLB over the z-position

FOV image:
1. run 'combine_files_and_apply_drift_correction': this takes the drift present in each frame, and updates the position of each localization.

2a. run 'compare_original_model': it produces the zoomed-in images of the FOV, and the colocalization. you can change the region that is zoomed in by changing the pixel limits (between 1 and 1024 pixels)

2b. run 'Run_combined_locations': Produces the full FOV images of the HeLa Cells.

Histograms (localizations):
1. run 'combine_files_and_apply_drift_correction', or use the file it returns
2. run'check_z_range_patches': it produces the histogram of the initial and second fit over z

Histograms (Outliers)
1. run 'check_outlier_positions': it takes the individual files as input

Z-displacement, Quiver plot, Chi-Squared: 
run 'analyze_fiterror_and_arrows'. Use one of the files, as the size of 'mu' and 'allspots' limits the size of data that can be accessed.

Table results:
run 'linklocs_3D_runlength_analysis': Note that you do either need to use a windows C compiler, or a Linux C compiler, to compile "cMultiFrameConnect_allcalcs3D".

All files or folders that are direct copies of 'vectorfit' have 'copy' added to their file or folder name. The 'vectorfit master' zip file and instructions are available on the 'vectorfit' repository. All rights are reserved to the Imphys Department

Rieger, B., Stallinga, S., Ligteringen, R., Droste, I., Martens, K., & Höppener, T. (2024). Vec-
torfit: Matlab implementation of the vector fitting method [Accessed: 2025-06-03].
https://gitlab.tudelft.nl/imphys/ci/vectorfit

The Dipimage library is used to visualize certain results. 
Luengo, C., & contributors. (2024). Diplib: Quantitative image analysis in c++, matlab and
python [Accessed: 2025-06-03]. https://github.com/DIPlib/diplib

The OTF fit was run in matlab:
Inc., T. M. (2024). Matlab version: 24.2.0.2833386 (r2024b). Natick, Massachusetts, United States.
https://www.mathworks.com

alongside certain toolboxes:
DIPimage, a scientific image analysis toolbox         Version 2.9                       License unknown

DIPlib, a scientific image analysis library           Version 2.9                       License unknown

Image Processing Toolbox                              Version 24.2        (R2024b)      License 329139 

Parallel Computing Toolbox                            Version 24.2        (R2024b)      License 329139 

Signal Processing Toolbox                             Version 24.2        (R2024b)      License 329139 

Statistics and Machine Learning Toolbox               Version 24.2        (R2024b)      License 329139 

Vehicle Network Toolbox                               Version 24.2        (R2024b)      License 329139 


P.S. if the OTF fit fails to properly run after the initial fit, remove the vectorfit master from your directory: certain files could be modified from previous versions. 

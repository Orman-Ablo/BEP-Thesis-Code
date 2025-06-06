# BEP-Thesis-Code
A collection of Code used to Produce the OTF fitting Routine used in the BEP thesis

In order to use The second fit, Follow the pipeline set for the 'vectorfit' Project:
https://gitlab.tudelft.nl/imphys/ci/vectorfit

Alternatively, the code for the 'Vectorfit' Project can be found in the 'vectorfit-master'.zip folder.


Once you have collected the initial fit of your localizations, follow the Pipeline for the OTF fit:

1. run the file 'combined_locations', which isolates the localizations of all ROIs
2. run 'PSF_data': it produces the 3-D stacks of the PSFs required for the 3-D OTF. you can change the upscaling factor to  produce finer/rougher sampling of the PSF.
3. run 'Plot_PSF_Results': this file produces the used 3-D OTF, and plots the shapes of the PSFs and OTFs. You can define here the desired size of the OTF
4. run 'localize_otfmode_remove_outliers': it is a modified copy of the 'localize'file used in the initial fit. It should be albe to run both the initial and second fit.

If you want to produce any of the plots used in the report, their pipeline is listed below:

CRLB:
1. run 'CRLB_calc': it combines the CRLB, and the positions
2. run 'plot_CRLB': Produces the figure of the mean CRLB over the z-position

FOV image:
1. run 'combined locations'
2a. run 'compare_original_model': it produces the zoomed-in images of the FOV, and the colocalization. you can change the region that is zoomed in by changing the pixel limits
2b. run 'Run_combined_locations': Produces the full FOV images of the HeLa Cells.

Histograms (localizations):
1. run'combined locations', or use the file it returns
2. run'check_z_range_patches': it produces the histogram of the initial and second fit over z

Histograms (Outliers)
1. run 'check_outlier_positions': it takes the individual files as input

Heatmap, Quiver plot, Chi-Squared: 
run 'analyze_fiterror_and_arrows'

Table results:
run 'linklocs_3D_runlength_analysis': Note that you do either need to use a windows C compiler, or a Linux C compiler.

All files or folders that are direct copies of 'vectorfit' have 'copy' added to their file or folder name. All rights are reserved to the Imphys Department

Rieger, B., Stallinga, S., Ligteringen, R., Droste, I., Martens, K., & Höppener, T. (2024). Vec-
torfit: Matlab implementation of the vector fitting method [Accessed: 2025-06-03].
https://gitlab.tudelft.nl/imphys/ci/vectorfit


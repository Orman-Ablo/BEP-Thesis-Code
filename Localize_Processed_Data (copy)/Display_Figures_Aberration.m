%Just to visualize the acquired data independently

clear all
close all

input = "C:/Users/Naam/Documents/BEP/Code/vectorfit-master/vectorfit-master/BEP Code Roman/Data_BEP_Roman/4_Aberration_Estimate\estimated_aberrations_rep_0001_init_0001_localization_segmentation_microtubules_3D_data0001.mat";
params = set_parameters_microtubules_3D;
loaded_data = load(input);
allspots = loaded_data.allspots;
roixy = loaded_data.roixy;
framelist = loaded_data.framelist;
theta_full = loaded_data.theta_




if params.show_plots
show_fitted_data
end
% This script is for fitting the local and global(aberration) parameters of
% synthetic single molecule data. This data can be generated using
% 'get_NAT_data.m'.
%
% INPUT
% allspots = the measured spots (segmented single molecule data)
% roixy = the locations of the measured spots in pixels.
% params = parameter file
% groundtruth = local and global parameters to compare to outcome
%
% OUTPUT
% -theta = estimated local and global variables
% -meritstore = the value of the cost function
% -mu = the modeled psfs
% -dmudtheta = derivative of mu to thetas
% through the iterations.
% -thetainit
%
% Autor: Isabel Droste, TU Delft, 2022, 2023, 2024

if params.simulated_data
    clearvars -except r params path_input_data path_input_data_full path_output_data path_output_data_full plot_data fit_without_aberrations nr_of_runs average_crlb_aberrations_total average_wrms_error_total max_crlb_aberrations_total max_wrms_error_total final_local_update path_gt_gammas system_aberration_error_all CRLB_all init nr_of_inits rep filename;
else
    clearvars -except r params allspots roixy framelist path_input_data path_input_data_full path_output_data path_output_data_full plot_data fit_without_aberrations nr_of_runs average_crlb_aberrations_total average_wrms_error_total max_crlb_aberrations_total max_wrms_error_total final_local_update init nr_of_inits rep filename;
end
close all 

%% Load data
loaded_data = load(path_input_data_full);
if params.simulated_data
    loaded_data = load(path_input_data_full);
    if isfield(loaded_data,'groundtruth')
        params.groundtruth_exists = true;
    else
        params.groundtruth_exists = false;
    end
    allspots = loaded_data.allspots;
    roixy = loaded_data.roixy;
    framelist = loaded_data.framelist;
    if params.groundtruth_exists
        groundtruth = loaded_data.groundtruth;
    end
    ID = loaded_data.ID;
else
    params.groundtruth_exists = false;
    fprintf('Load input data\n')    
    fprintf('Finished loading input data\n')
    allspots = loaded_data.allspots;
    roixy = loaded_data.roixy;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(roixy,2) == 2
        roixy = permute(roixy,[2 1]);
    end
    framelist = loaded_data.localizations_with_outliers(:,5);
    ID = loaded_data.localizations_with_outliers(:,1);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% select spots
params.Ncfg_total = size(allspots,4);
if params.select_spots
    if init==1
        fprintf('Select %i spots for fitting\n',params.Ncfg);
        %params.selected_indices = sort(randperm(params.Ncfg_total,params.Ncfg)); %random selection
        params.selected_indices = select_spots(roixy,params); % selection on a grid for uniform distribution
    else
        % Use spots from rep 1
    end
else
    params.Ncfg = params.Ncfg_total;
    params.selected_indices = 1:params.Ncfg;
    params.Ncfg_global = params.Ncfg;
end


%%%%%%
% Determine the random indices for each global update, so we can save and
% repeat them
if params.select_spots
    if init==1
        Ncfg_global = params.Ncfg_global;
        Ncfg = params.Ncfg;
        params.selected_indices_global_all = zeros(Ncfg_global,params.max_total_iterations);
        for i = 1:params.max_total_iterations
            selected_indices_global_all(:,i) = randperm(Ncfg,Ncfg_global);
        end
        params.selected_indices_global_all = selected_indices_global_all;
    else
        % Use spots from rep 1
    end
end

% remove unselected spots
allspots_original = allspots;
roixy_original = roixy;
indices = params.selected_indices;
allspots = allspots(:,:,:,indices);
roixy = roixy(:,indices);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framelist = framelist(indices);
ID = ID(indices);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.groundtruth_exists
    groundtruth_original = groundtruth;
    groundtruth.local = groundtruth.local(:,indices);
end

% % Show selected spots
% figure
% scatter(roixy(1,:),roixy(2,:),20,'filled', 'MarkerFaceAlpha', 0.5)
% axis square
% pause

%% Fit from first initial value

[thetainit,theta_temp,thetastore_temp,meritstore_temp,mu_temp,dmudtheta_temp,Niters_temp,outliers,localizations,localizations_with_outliers,monitorstore_total] = localization_nat(allspots,roixy,framelist,ID,params);

nr_of_initial_values = 1;
theta_full.local = zeros([size(theta_temp.local),nr_of_initial_values]);
theta_full.global = zeros([size(theta_temp.global),nr_of_initial_values]);
thetastore_full.local = zeros([size(thetastore_temp.local),nr_of_initial_values]);
thetastore_full.global = zeros([size(thetastore_temp.global),nr_of_initial_values]);
localizations_full = zeros([size(localizations),nr_of_initial_values]);
localizations_with_outliers_full = zeros([size(localizations_with_outliers),nr_of_initial_values]);
meritstore_full.local = zeros([size(meritstore_temp.local),nr_of_initial_values]);
meritstore_full.global = zeros([size(meritstore_temp.global),nr_of_initial_values]);
mu_full = zeros([size(mu_temp),nr_of_initial_values]);
dmudtheta_full.local = zeros([size(dmudtheta_temp.local),nr_of_initial_values]);
dmudtheta_full.global = zeros([size(dmudtheta_temp.global),nr_of_initial_values]);
Niters_full.local = zeros([size(Niters_temp.local),nr_of_initial_values]);
Niters_full.global = zeros([size(Niters_temp.global),nr_of_initial_values]);
outliers_full = cell(nr_of_initial_values,1);

theta_full.local(:,:,1) = theta_temp.local;
theta_full.global(:,:,1) = theta_temp.global;
thetastore_full.local(:,:,:,1) = thetastore_temp.local;
thetastore_full.global(:,:,1) = thetastore_temp.global;
localizations_full(:,:,1) = localizations;
localizations_with_outliers_full(:,:,1) = localizations_with_outliers;
meritstore_full.local(:,:,1) = meritstore_temp.local;
meritstore_full.global(:,:,1) = meritstore_temp.global;
mu_full(:,:,:,:,1) = mu_temp;
dmudtheta_full.local(:,:,:,:,:,1) = dmudtheta_temp.local;
dmudtheta_full.global(:,:,:,:,:,1) = dmudtheta_temp.global;
Niters_full.local(:,:,1) = Niters_temp.local;
Niters_full.global(:,:,1) = Niters_temp.global;
outliers_full{1} = outliers;

merit_last_full = zeros(nr_of_initial_values,1);
merit_max = meritstore_temp.global(end);
merit_last_full(1) = meritstore_temp.global(end);
max_merit_index = 1;

%% Fit again from other initual values
if params.multiple_initial_gammas
    for s=2:nr_of_initial_values
        fprintf('Start fitting initial value %i/%i\n',s,nr_of_initial_values);
        
        thetainit.global = thetainit_all.global(:,s);
        thetainit.local = initialvalues_phasor(allspots,roixy,thetainit.global,params);
        thetainit_all.local(:,:,s) = thetainit.local;
        
        [theta_temp,thetastore_temp,meritstore_temp,mu_temp,dmudtheta_temp,Niters_temp,outliers,localizations,localizations_with_outliers] = localization_nat(allspots,roixy,thetainit,framelist,ID,params);

        theta_full.local(:,:,s) = theta_temp.local;
        theta_full.global(:,:,s) = theta_temp.global;
        thetastore_full.local(:,:,:,s) = thetastore_temp.local;
        thetastore_full.global(:,:,s) = thetastore_temp.global;
        localizations_full(:,:,s) = localizations;
        localizations_with_outliers_full(:,:,s) = localizations_with_outliers;
        meritstore_full.local(:,:,s) = meritstore_temp.local;
        meritstore_full.global(:,:,s) = meritstore_temp.local;
        mu_full(:,:,:,:,s) = mu_temp;
        dmudtheta_full.local(:,:,:,:,:,s) = dmudtheta_temp.local;
        dmudtheta_full.global(:,:,:,:,:,s) = dmudtheta_temp.global;
        Niters_full.local(:,:,s) = Niters_temp.local;
        Niters_full.global(:,:,s) = Niters_temp.global;
        outliers_full{s} = outliers;
    
        merit_last_full(s) = meritstore_temp.global(end);

        if merit_last_full(s) > merit_max
            max_merit_index = s;
            merit_max = merit_last_full(s);
        end
    
    end
end

disp('Estimation aberrations completed')
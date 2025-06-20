% Try out: select winning solution based on other spots that it has not
% seen before.

close all
clear all

path_input_data = '/home/idroste/Desktop/TUDelft/Data/Lidke/data_analysis/DATA1/Data_wo_Astigmatism/segmentation/segmentation_thr60_Data0001.mat';
path_fitted_data = '/home/idroste/Desktop/TUDelft/Data/Lidke/data_analysis/DATA1/Data_wo_Astigmatism/fitted_aberrations/rand_init_fixed_merit/0*';

% Load full dataset
path_input_data_full = path_input_data;
load(path_input_data_full);

params = set_parameters_Lidke_2D;
params.Ncfg_total = size(allspots,4);

% Select reference spots
N_ref = 1000;
selected_indices = sort(randperm(params.Ncfg_total,N_ref));

allspots_original = allspots;
roixy_original = roixy;
allspots = allspots(:,:,:,selected_indices);
roixy = roixy(:,selected_indices);
framelist = framelist(:,selected_indices);
ID = ID(:,selected_indices);

% Test solutions one by one
all_solutions = dir(path_fitted_data);

merit_all = zeros(size(all_solutions,1),1);
merit_all_no_outliers = zeros(size(all_solutions,1),1);
merit_all_outliers = zeros(size(all_solutions,1),1);

for i=1:size(all_solutions,1)
    path_full = fullfile(all_solutions(i).folder,all_solutions(i).name);
    load(path_full,'theta_full')
    theta_global_temp = theta_full.global;

    % do local update
    numparams = params.numparams;
    numgammas = params.numgammas;
    Ncfg = N_ref;
    params.Ncfg = Ncfg;
    min_total_iterations = params.min_total_iterations;
    max_total_iterations = params.max_total_iterations;

    thetastore_local = zeros(numparams,Ncfg,0);
    thetastore_global = zeros(numgammas,0);
    
    meritstore_local = zeros(Ncfg,0);
    meritstore_global = [];
    monitorstore_total = [];
    
    alambda_local = ones(Ncfg,1)*params.alambda0_local;
    alambda_global = params.alambda0_global;
    
    Niters_local = zeros(Ncfg,0);
    Niters_global = [];
    
    params.perform_final_iterations = false;
    flip_z_net = false;
    thetainit.global = theta_global_temp;
    if params.cpp
        fitmodel_matlab = params.fitmodel;
        if strcmp(fitmodel_matlab,'xy-gamma')
            params.fitmodel = 'xy-constaberrations';
        elseif strcmp(fitmodel_matlab,'xyz-gamma')
            params.fitmodel = 'xyz-constaberrations';
        end
        fitter = vectorfitter(allspots,roixy, params,params.cpp_fitmode);
        fitter.recompute_zernike_coefficients(thetainit.global, params.aberrations);
        thetainit.local = fitter.theta_init();
        fitter.recompute_zernike_coefficients(thetainit.global, params.aberrations);
        params.fitmodel = fitmodel_matlab;
    else
        thetainit.local = initialvalues_phasor(allspots,roixy,thetainit.global,params);
        % not sure yet of the effect of this.
        %thetainit.local = thetainit.local+0.1*(rand(params.numparams,params.Ncfg)-0.5)*params.pixelsize;
    end
    
    theta_local = thetainit.local;
    theta_global = thetainit.global
    iiter_total = 1;

    if params.cpp
        [fitter,theta_local,thetastore_local,localizations,localizations_with_outliers,meritstore_local,mu,dmudtheta_local,Niters_local,outliers] = ...
        local_update_cpp(fitter,thetastore_local,theta_global,meritstore_local,allspots,roixy,iiter_total,Niters_local,flip_z_net,framelist,ID,params); 
    else
        [theta_local,thetastore_local,localizations,localizations_with_outliers,meritstore_local,alambda_local,mu,dmudtheta_local,Niters_local,outliers] = ...
        local_update(theta_local,thetastore_local,theta_global,meritstore_local,alambda_local,allspots,roixy,iiter_total,Niters_local,flip_z_net,framelist,ID,params);    
    end

    %for jcfg = 1:Ncfg
    %    merit_all(i) = merit_all(i) + likelihood(params,allspots(:,:,:,jcfg),mu(:,:,:,jcfg),dmudtheta_local(:,:,:,:,jcfg),params.varfit);
    %end

    no_outliers = setdiff(1:Ncfg,outliers);
    merit_all(i) = mean(meritstore_local);
    merit_all_no_outliers(i) = mean(meritstore_local(no_outliers));
    merit_all_outliers(i) = mean(meritstore_local(outliers));


end

%% Plot likelihoods

figure
x = merit_all;
bar(x)
ylim([min(x)-0.1*(max(x) - min(x)) max(x)+0.1*(max(x) - min(x))])

%%
figure
x = merit_all_no_outliers;
bar(x)
ylim([min(x)-0.1*(max(x) - min(x)) max(x)+0.1*(max(x) - min(x))])

%%
figure
x = merit_all_outliers;
bar(x)
ylim([min(x)-0.1*(max(x) - min(x)) max(x)+0.1*(max(x) - min(x))])

%%
mean(theta_local(4,outliers))
mean(theta_local(4,no_outliers))
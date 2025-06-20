% Fit all ROIs with (or without) the estimated aberrations.
% The input data should be segmented
% You can choose between GPU, CPU or Matlab fitting

clear all
close all

input_path = 'C:/Users/Naam/Documents/BEP/Code/vectorfit-master/vectorfit-master/data/1_segmentation/segmentation_All_framesmicrotubules_3D_data0001.mat';

path_aberrations = 'C:/Users/Naam/Documents/BEP/Code/vectorfit-master/vectorfit-master/data/2_estimate_aberrations/estimated_aberrations_All_framesrep_0001_init_0001_segmentation_microtubules_3D_data0001.mat';

output_path = 'C:/Users/Naam/Documents/BEP/Code/vectorfit-master/vectorfit-master/data/3_localize/';

params = set_parameters_microtubules_3D;

params.use_fitted_aberrations = true; % If false, aberrations are set to 0, except for constant astigmatism for 3D data
params.path_aberrations = path_aberrations;

params.fitmodel = 'xyz';
spots_per_batch = 2e4;

if params.cpp
    if strcmp(params.fitmodel,'xy')
        params.fitmodel = 'xy_constaberrations';
    elseif strcmp(params.fitmodel,'xyz')
        params.fitmodel = 'xyz_constaberrations';
    end
end

all_input_files = dir(input_path);
nr_of_files = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_files);

tic_total = tic;

for file_i=1:nr_of_files
    
    % Load data
    filename = all_input_files(file_i).name; 
    foldername = all_input_files(file_i).folder;
    path_input_data_full = fullfile(foldername, filename);    
    fprintf('Input segmented data = %s\n',path_input_data_full);

    loaded_data = load(path_input_data_full);
    allspots = loaded_data.allspots;
    roixy = loaded_data.roixy;
    framelist = loaded_data.framelist;
    ID = loaded_data.ID;

    params.select_spots = false;
    Ncfg_total = size(allspots,4);
    params.Ncfg_total = Ncfg_total;
    params.Ncfg_global = Ncfg_total;
    params.Ncfg = Ncfg_total;
    
    fprintf('\nPerform localization with %i ROIs\n',Ncfg_total);
    
    max_local_iterations = params.max_local_iterations;
    numparams = params.numparams;

    if params.use_fitted_aberrations
        fitted_aberrations_data = load(path_aberrations,'theta_full','max_merit_index');
        theta_global = fitted_aberrations_data.theta_full.global;
    else    
        theta_global = zeros(params.numgammas,1);  
        theta_global(6) = 66; % Microtubules 3D average value
    end

    Nbatch = ceil(Ncfg_total/spots_per_batch);

    theta_init = zeros(params.numparams,Ncfg_total);
    outliers = [];
    mu = [];
    merit = zeros(Ncfg_total,1);
    meritoffset = zeros(Ncfg_total,1);
    Niters = [];

    % ID x y z framenumber CRLBx CRLBy CRLBz Nph bg
    localizations_with_outliers = zeros(Ncfg_total,10);
    localizations_with_outliers(:,1) = ID;
    localizations_with_outliers(:,5) = framelist;

    % Precalculation of OTFs
    if params.cpp && params.FlagOTF
        fitter_all = vectorfitter(allspots,roixy,params,params.cpp_fitmode);
        [otfs,spot_otf_indices] = fitter_all.update_otfs(theta_global,params.aberrations);
        delete(fitter_all)
    end
    
    % Fit in batches
    tic_file = tic;
    for batch = 1:Nbatch
        batch
        start = (batch-1)*spots_per_batch;
        stop = min(start+spots_per_batch,Ncfg_total);
        Nspots = stop - start;

        spots_batch = allspots(:,:,:,start+1:stop);
        roixy_batch = roixy(:,start+1:stop);

        if params.cpp
            fitter = vectorfitter(spots_batch,roixy_batch, params,params.cpp_fitmode);
            if params.FlagOTF
                spot_otf_indices_batch = spot_otf_indices(start+1:stop);
                fitter.set_otfs(otfs,spot_otf_indices_batch);
            else
                cpp_zernike_coeffs = fitter.update_zernike_coefficients(theta_global, params.aberrations);
            end
            theta_init_batch = fitter.theta_init();
            [theta_update_batch, Niters_batch] = fitter.local_update();
            dmudtheta_batch = fitter.dmudtheta;
            outliers_batch = fitter.outliers();
            crlb = fitter.fisher_crlb();
            mu_batch = fitter.mu();
            mu = cat(4,mu,mu_batch); 
            delete(fitter)
        else
            Ncfg = size(roixy_batch,2)
            params.Ncfg = Ncfg;
            params.Ncfg_total = Ncfg;
            thetastore_local = zeros(numparams,Ncfg,0);        
            meritstore_local = zeros(Ncfg,0);   
            alambda_local = ones(Ncfg,1)*params.alambda0_local;         
            Niters_batch = zeros(Ncfg,0);   
            iiter_total = 1;
            params.perform_final_iterations = false;
            flip_z_net = false;
            theta_init_batch = initialvalues_phasor(spots_batch,roixy_batch,theta_global,params);
            [theta_update_batch,thetastore_local,localizations_temp,localizations_with_outliers_temp,meritstore_local,alambda_local,mu_batch,dmudtheta_batch,Niters_batch,outliers_batch] = ...
            local_update(theta_init_batch,thetastore_local,theta_global,meritstore_local,alambda_local,spots_batch,roixy_batch,iiter_total,Niters_batch,flip_z_net,framelist(start+1:stop),ID(start+1:stop),params); 
            crlb = localizations_with_outliers_temp(:,6:8)';
            outliers_batch = outliers_batch';
        end

        % Calculate the merit
        Ncfg = size(theta_update_batch,2);
        for jcfg = 1:Ncfg
            merit(start+jcfg) = likelihood(params,allspots(:,:,:,start+jcfg),mu_batch(:,:,:,jcfg),dmudtheta_batch(:,:,:,:,jcfg),params.varfit);
        end
        meritoffset(start+1:stop) = meritoffsetcalc(allspots(:,:,:,start+1:stop),params.varfit);

        no_outliers_temp = setdiff(1:Nspots,outliers_batch);
        outliers = cat(1,outliers,outliers_batch+start);
        theta_init(:,start+1:stop) = theta_init_batch;

        [~,xy_um] = get_fov_coordinates(roixy_batch,theta_update_batch(1,:),theta_update_batch(2,:),params);
        
        localizations_with_outliers(start+1:stop,2:3) = 1e3*xy_um';
        if contains(params.fitmodel,'xyz')
            localizations_with_outliers(start+1:stop,4) = theta_update_batch(3,:);
            localizations_with_outliers(start+1:stop,9:10) = theta_update_batch(4:5,:)';
            localizations_with_outliers(start+1:stop,6:8) = crlb(1:3,:)';
        else %xy
            localizations_with_outliers(start+1:stop,4) = zeros(Nspots,1);
            localizations_with_outliers(start+1:stop,9:10) = theta_update_batch(3:4,:)';
            localizations_with_outliers(start+1:stop,6:7) = crlb(1:2,:)';
            localizations_with_outliers(start+1:stop,8) = zeros(Nspots,1);
        end
        
        theta_update(:,start+1:stop) = theta_update_batch;
        theta_init(:,start+1:stop) = theta_init_batch;
        mu(:,:,:,start+1:stop) = mu_batch;
        Niters(start+1:stop) = Niters_batch;

    end

    merit_corrected = merit+meritoffset;
    no_outliers = setdiff(1:Ncfg_total,outliers);
    localizations = localizations_with_outliers(no_outliers,:);

    t_file = toc(tic_file)
    locs_per_second = Ncfg_total/t_file

    filename = all_input_files(file_i).name; 
    path_output_data_full = strcat(output_path,'localization_',filename);
    save(path_output_data_full,'theta_update','theta_global','localizations','localizations_with_outliers','theta_init','outliers','roixy','merit','merit_corrected','params','allspots','mu','-v7.3')
    %save(path_output_data_full,'theta_update','theta_global','localizations','localizations_with_outliers','theta_init','outliers','roixy','merit','merit_corrected','params','-v7.3')
    
    fprintf('Localization with %i ROIs finished\n',Ncfg_total)

end

t_total = toc(tic_total)

%% Render final image
if params.show_plots
    x1 = localizations(:,2);
    y1 = localizations(:,3);
    z1  = localizations(:,4);
    
    locations = [x1 y1 z1];
   
    ncolors = 19;
    colorrange = [-500,400];
    
    xlim = [0.0 1.0];
    ylim = [0.0 1.0];
    pixelsize = 50

    sigma = 2.0;
    clim = 0.95; % the brightest (1-clim) fraction of pixels gets clipped 
    render_3D_func(locations,ncolors,colorrange,xlim,ylim,pixelsize,sigma,clim,params);
end


%% Combine localizations from multiple files into one file
% 
% input_folder_combine = strcat(output_path,'localization*');
% output_path_final = strcat(output_path,'localization_combined.mat');
% 
% all_input_files = dir(input_folder_combine);
% nr_of_files = size(all_input_files,1);
% fprintf('Found %i files\n',nr_of_files);
% 
% localizations = [];
% %localizations_with_outliers = [];
% %theta_init = [];
% %theta = [];
% %outliers = [];
% roixy = [];
% %allspots = [];
% %mu = [];
% Ncfg_total = 0;
% 
% for file_i=1:nr_of_files
% %for file_i=1:7
% 
%     % Load data
%     filename = all_input_files(file_i).name; 
%     foldername = all_input_files(file_i).folder;
%     path_input_data_full = fullfile(foldername, filename);    
%     fprintf('Input data = %s\n',path_input_data_full);
%     loaded_data = load(path_input_data_full);
% 
%     localizations_temp = loaded_data.localizations;
%     %localizations_with_outliers_temp = loaded_data.localizations_with_outliers;
%     %theta_init_temp = loaded_data.theta_init;
%     %outliers_temp = loaded_data.outliers;
%     roixy_temp = loaded_data.roixy;
%     %allspots_temp = loaded_data.allspots;
%     %mu_temp = loaded_data.mu;
%     %theta_temp = loaded_data.theta_update;
%     params = loaded_data.params;
% 
%     localizations = cat(1,localizations,localizations_temp);
%     %localizations_with_outliers = cat(1,localizations_with_outliers,localizations_with_outliers_temp);
%     %theta_init = cat(2,theta_init,theta_init_temp);
%     %outliers_temp_adjusted = outliers_temp + Ncfg_total;
%     %outliers = cat(1,outliers,outliers_temp_adjusted);
%     %Ncfg_total = Ncfg_total + size(roixy_temp,2);
%     roixy = cat(2,roixy,roixy_temp);
%     %allspots = cat(4,allspots,allspots_temp);
%     %mu = cat(4,allspots,allspots_temp);
%     %theta = cat(2,theta,theta_temp);
% end
% 
% save(output_path_final,'localizations','roixy','params')



%Notes:
%Why is spots_per batch 2e4?
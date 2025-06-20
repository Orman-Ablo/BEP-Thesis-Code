% This script plots multiple bead spots and if desired the corresponding 
% fitted psfs.
% For each bead, a z-profile can be plotted. That is a plot of the sum of
% the psf for each z-value.

% Input
% - Spots: synthetically genereated or segemented bead data (.m file).
% - Localization results including the parameter 'mu' that contains the
% fitted (N * PSF + b).

% Output
% - 3D plot of the beads and fitted PSF
% - z-profile that comapares the total intensity per z-slice of the data
% and the fit.

addpath('C:\Users\idroste\Desktop\TUDelft\Code\diplib\share\DIPimage');
setenv('PATH',['C:\Users\idroste\Desktop\TUDelft\Code\diplib\bin',';',getenv('PATH')]);

data_in_workspace = false;
include_fitted_data = false;
include_outfocus = true;
scale_images = false; % Scale images such that max of data and fits is the same, for better visualization.

if ~data_in_workspace
    clearvars -except data_in_workspace include_outfocus include_fitted_data scale_images;
end

%% Input data

%source_path_spots = 'N:\tnw\IST\QI\users\idroste\Data\data_bead_fitting_22_02_15\processed_archive\50iterations\zernikes12\segmentation_old_region_3.mat';
%source_path_localization = 'N:\tnw\IST\QI\users\idroste\Data\data_bead_fitting_22_02_15\processed_archive\50iterations\zernikes12\localization_old_region_3.mat';

% Synthetic data
%source_path_spots = 'N:\tnw\IST\QI\users\idroste\Data\synthetic_data\small_beads\bead_data.mat';
%source_path_localization = 'N:\tnw\IST\QI\users\idroste\Data\synthetic_data\small_beads\results\localization_bead_data.mat';
%source_path_spots = 'N:\tnw\IST\QI\users\idroste\Data\synthetic_data\beads_aberrations4\bead_data.mat';
%source_path_localization = 'N:\tnw\IST\QI\users\idroste\Data\synthetic_data\beads_aberrations4\results_centroid\localization_bead_data.mat';


% Real data
source_path_spots = 'N:\tnw\IST\QI\users\idroste\Data\data_bead_fitting_22_02_15\processed\zernike12_50it\segmentation_old_region_5';
source_path_localization = 'N:\tnw\IST\QI\users\idroste\Data\data_bead_fitting_22_02_15\processed\zernike12_50it\localization_old_region_5';

%source_path_spots = 'C:\Users\idroste\Desktop\TUDelft\Code\vecfitcpu_vortex\read\segmentation_old_region_5.mat';
%source_path_localization = 'C:\Users\idroste\Desktop\TUDelft\Code\vecfitcpu_vortex\read\localization_old_region_5.mat'

%destination_path = 'C:\Users\idroste\Desktop\TUDelft\Data\data_bead_fitting_22_02_15\processed\';

% The bead of which the data will be saved.
beads = 76:76;


%% Select the data

% Open the input data.
if ~data_in_workspace
    spots_data = load(source_path_spots);
    fprintf('Spots data loaded\n');
    localization_results = load(source_path_localization);
    fprintf('Localization data loaded\n');
end

% Get the measured data for the selected beads.
all_measured_beads = spots_data.allspots(:,:,:,:);
params_measured = spots_data.params; 
%all_measured_beads = (all_measured_beads-params_measured.offset)/params_measured.gain;
%all_measured_beads(all_measured_beads<=0) = 1e-3; 
if ~include_outfocus
    outfocus = localization_results.outfocus;
    all_measured_beads(:,:,:,outfocus) = [];
end

% Create image with beads next to each other.
measured_beads_img = all_measured_beads(:,:,:,beads(1));
for i = 2:length(beads)
    measured_beads_img = cat(2,measured_beads_img,all_measured_beads(:,:,:,beads(i)));
end

% Create images of the fitted PSFs
if include_fitted_data
    params = localization_results.params;  
    theta_end = localization_results.thetastore_original(:,:,end); % params, beads, iters
    [fitted_psfs,~] = poissonrate(params,theta_end(:,beads(1)));
    for i = 2:length(beads)
        disp(i);
        [psf,~] = poissonrate(params,theta_end(:,beads(i)));
        fitted_psfs = cat(4,fitted_psfs,psf);
    end

    % Create image of fitted PSFs
    theta_intermediate = localization_results.thetastore_original(:,1,end);
    [fitted_psfs, ~] = poissonrate(params,theta_intermediate);
    for i = 2:length(beads)
        disp(i);
        if include_outfocus
            theta_intermediate = localization_results.thetastore_original(:,i,end);
        else
            theta_intermediate = localization_results.thetastore(:,i,end);
        end
        [intermediate_psf,~] = poissonrate(params,theta_intermediate);
        fitted_psfs = cat(2,fitted_psfs,intermediate_psf);
    end
    
    % Scale images.
    if scale_images
        max_beads = max(measured_beads_img, [], 'all');
        max_fits = max(fitted_psfs, [], 'all');
        fitted_psfs = (max_beads/max_fits)*fitted_psfs;
    end

    all_data = cat(1,measured_beads_img,fitted_psfs);
    
else
    all_data = measured_beads_img;
end % end if include_fitted_data

% Show measured and fitted beads.
dipshow(measured_beads_img,'lin');
diptruesize(500)
colormap parula

dipshow(fitted_psfs,'lin');
diptruesize(500)
colormap parula

if false 
    % Create z-profiles
    z_range = params.zrange;
    Mz = params.Mz;
    z = 1e-3*linspace(z_range(1),z_range(2),Mz);
    
    total_sum_beads = zeros(length(beads),1);
    total_sum_psfs = zeros(length(beads),1);
    
    for i = 1:length(beads)
        selected_bead = all_measured_beads(:,:,:,beads(i));
        selected_psf = fitted_psfs(:,:,:,i);
    
        xy_sum_bead = squeeze(sum(selected_bead,[1 2]));
        xy_sum_psf = squeeze(sum(selected_psf,[1 2]));
    
        total_sum_beads(i) = squeeze(sum(xy_sum_bead,"all"));
        total_sum_psfs(i) = squeeze(sum(xy_sum_psf,"all"));
    
        figure
        plot(z,xy_sum_bead,'LineWidth',1.5);
        hold on
        plot(z,xy_sum_psf,'LineWidth',1.5);
        legend('Measured bead','Fitted PSF');
        xlabel('z (\mu m)');
        ylabel('sum intensity');
    end
end


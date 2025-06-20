% Plot merit landscape of experimental data

close all
clear all

%path_input_data = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Yutong_Nanorulers\data_analysis\CMOS\nanorulers\segmentation\segmentation_NR80R_1201_Oil100x_Col175_647P300-100_100ms_2000f_FL_Pos2.mat';
%path_input_data = '/home/idroste/Desktop/TUDelft/Code/temp/smlm-vector-localization/read\simulationdata\input_data\input_spots_nanorulerbeads_without_tilt_1_.mat';
path_input_data = 'C:\Users\idroste\Desktop\TUDelft\Code\smlm-branches\smlm-vector-localization\read\simulationdata\input_data\input_spots_astig_coma_1_.mat';

%path_fitted_data = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Yutong_Nanorulers\data_analysis\CMOS\nanorulers\fitted_aberrations\fitted_aberrations_NR80R_1201_Oil100x_Col175_647P300-100_100ms_2000f_FL_Pos2_018.mat';
%path_fitted_data = '/home/idroste/Desktop/TUDelft/Code/temp/smlm-vector-localization/read\simulationdata\results\results_nanorulerbeads_without_tilt_1_.mat';
path_fitted_data = 'C:\Users\idroste\Desktop\TUDelft\Code\smlm-branches\smlm-vector-localization\read\simulationdata\results\results_nanorulerbeads_astig_coma_1_.mat';

fitted_data = load(path_fitted_data);
theta_local = fitted_data.theta_full.local;
theta_global = fitted_data.theta_full.global;
selected_indices = fitted_data.params.selected_indices;
groundtruth = fitted_data.groundtruth;
%theta_global = groundtruth.global;
%theta_local = groundtruth.local;

input_data = load(path_input_data);
allspots = input_data.allspots(:,:,:,selected_indices);
roixy = input_data.roixy(:,selected_indices);
params = input_data.params;
params.nr_of_parallel_workers = 8;

params.selected_indices_global = 1:200;

%%
numparams = params.numparams;
numgammas = params.numgammas;
Ncfg = params.Ncfg;

thetastore_local = zeros(numparams,Ncfg,0);
thetastore_global = zeros(numgammas,0);

meritstore_local = zeros(Ncfg,0);
meritstore_global = [];
monitorstore_total = [];

alambda_local = ones(Ncfg,1)*params.alambda0_local;
alambda_global = params.alambda0_global;

Niters_local = zeros(Ncfg,0);
Niters_global = [];

iiter_total = 1;
framelist = ones(Ncfg,1);
ID = 1:Ncfg;

flip_z_net = false;

start_x = -20;
stop_x = 5;
steps_x = 15;
x = linspace(start_x,stop_x,steps_x);

start_y = -25;
stop_y = 5;
steps_y = 15;
y = linspace(start_y,stop_y,steps_y);

[X,Y] = meshgrid(x,y);
merit_matrix = zeros(steps_y,steps_x);

%theta_local = initialvalues_phasor(allspots,roixy,theta_global,params);
%[theta_local,thetastore_local,localizations,localizations_with_outliers,meritstore_local,alambda_local,mu_allspots,dmudtheta_allspots,Niters_local,outliers] = ...
% local_update(theta_local,thetastore_local,theta_global,meritstore_local,alambda_local,allspots,roixy,iiter_total,Niters_local,flip_z_net,framelist,ID,params);


for xi=1:steps_x
    for yi=1:steps_y

        fprintf('xi,yi = %i,%i\n',xi,yi);

        % Set gammas
        theta_global(9) = x(xi);
        theta_global(12) = y(yi);

        %theta_local = initialvalues_phasor(allspots,roixy,theta_global,params);
        %[theta_local,thetastore_local,localizations,localizations_with_outliers,meritstore_local,alambda_local,mu_allspots,dmudtheta_allspots,Niters_local,outliers] = ...
        % local_update(theta_local,thetastore_local,theta_global,meritstore_local,alambda_local,allspots,roixy,iiter_total,Niters_local,flip_z_net,framelist,ID,params);

        [merit_total,merit_allspots,grad_total,Hessian_total,dmudtheta_allspots] = ...
        get_merit_nat(theta_local,theta_global,allspots,roixy,params);
        
        merit_matrix(yi,xi) = merit_total;

    end
end
%%
%merit_matrix_local_gt = merit_matrix;
%merit_matrix_crop_local_gt = merit_matrix_crop;
figure
merit_matrix_crop_local_gt = max(merit_matrix_crop,1109400);
surf(X,Y,merit_matrix_crop_fitted,'FaceAlpha',0.5)
hold on
surf(X,Y,merit_matrix_crop_local_gt,'FaceAlpha',0.5)
%ylim([-25 0])
zlim([1109400 1110200])
xlabel('gamma 9')
ylabel('gamma 12')
axis square
colorbar
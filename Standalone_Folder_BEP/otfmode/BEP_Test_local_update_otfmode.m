%% Testing the parameters needed for, and effectiveness of, the poissonrate using 'get_psfs_derivatives_otfmode' to acquire the PSF and its derivatives
clear all
close all
input_path ="C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_Edges.mat";
all_input_files = dir(input_path);
nr_of_datasets = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_datasets);

Example_file_spots = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Extra input data\transfer_3055388_files_4e4cf1e2\localization_segmentation_roi17_thr15p5_Data0001.mat";

% %% set up parallel pool
% nr_of_parallel_workers = 8;
% if isempty(gcp('nocreate'))
%     parpool('Processes', nr_of_parallel_workers);
% end




%% Parameters: Use one of the parameters to initialize some values (e.g. file 1)
%% Input (Averaged PSFs and OTFs)
dataset_name = all_input_files.name; 
foldername = all_input_files(1).folder;
input_filename = fullfile(foldername, dataset_name);
loaded_input = load(input_filename);
PSFmodelmatrix= loaded_input.PSFmodelmatrix;
load(input_path,"Xpatch","Ypatch",'Zpatch')

chisquarematrix= loaded_input.chisquarematrix;
numlocmatrix= loaded_input.numlocmatrix;
PSFdatamatrix= loaded_input.PSFdatamatrix;
OTF_3d = loaded_input.OTF_3d;
OTF_3d_abs = loaded_input.OTF_3d_abs;

all_spots_info = dir(Example_file_spots);
Spots_dataset = all_spots_info(1).name;
Spots_folder = all_spots_info(1).folder;
Spots_filename = fullfile(Spots_folder, Spots_dataset);
load_spots = load(Spots_filename);
localizations = load_spots.localizations;

framelist = localizations(:,5);
ID = localizations(:,1);



theta_update = load_spots.theta_update;
outliers = load_spots.outliers;
allspots = load_spots.allspots;
roixy = load_spots.roixy;
params = load_spots.params;
Npatch = size(PSFmodelmatrix,5);%params.Npatch;
params.Npatch =Npatch;
theta_global = load_spots.theta_global;
numparams = params.numparams;


locindex_file = 1:size(allspots,4);
masklocs = ~ismember(locindex_file,outliers);
%% Makes Sure all the used variables are properly Seperated!
Spots_fit = allspots(:,:,:,masklocs);
roixy_fit = roixy(:,masklocs);

theta_fit = theta_update(:,masklocs);


z_fit = [min(Zpatch) max(Zpatch)];
sep = find((localizations(:,4)>= z_fit(1)) & (localizations(:,4)<= z_fit(2)));
ID = ID(sep);
framelist = framelist(sep);
Spots_fit = Spots_fit(:,:,:,sep);
roixy_fit = roixy_fit(:,sep);
theta_fit = theta_fit(:,sep);

theta_init = zeros(params.numparams,Ncfg_total);
outliers = [];
mu = [];
merit = zeros(Ncfg_total,1);
meritoffset = zeros(Ncfg_total,1);
Niters = [];




% make sure all these parameters are present!
zemit = 0

params.Notf = params.Notfx;
params.Mroix = params.Mx;
params.Mroiy = params.My;
params.roisamplingdistance = params.pixelsize; % sampling distance in ROI
params.xroirange = params.roisamplingdistance*params.Mroix/2; % 1/2 of size PSF space in x
params.yroirange = params.roisamplingdistance*params.Mroiy/2; % 1/2 of size PSF space in y
params.fitmodel = 'xyz';
params.zstage = zemit;%% ask whther this makes sense
params.compders =1;
params.K=1;
params.zspread = [-1000,1000];



framelist = localizations_with_outliers(:,5);
spots_per_batch = 2e4;
ID = load_spots.localizations_with_outliers(:,1);
framelist = load_spots.localizations_with_outliers(:,5);
Ncfg_total = 2*spots_per_batch

% ID x y z framenumber CRLBx CRLBy CRLBz Nph bg
localizations_with_outliers = zeros(Ncfg_total,10);
localizations_with_outliers(:,1) = ID;
localizations_with_outliers(:,5) = framelist;
Nbatch = ceil(Ncfg_total/spots_per_batch);

thetatry = theta_update;
for batch = 1:Nbatch
    batch
    start = (batch-1)*spots_per_batch;
    stop = min(start+spots_per_batch,Ncfg_total);
    Nspots = stop - start;
    spots_batch = allspots(:,:,:,start+1:stop);
    roixy_batch = roixy(:,start+1:stop);
    theta_batch = theta_fit(:,start+1:stop);

    Ncfg = size(roixy_batch,2);
    params.Ncfg = Ncfg;
    params.Ncfg_total = Ncfg;
    thetastore_local = zeros(numparams,Ncfg,0);        
    meritstore_local = zeros(Ncfg,0);   
    alambda_local = ones(Ncfg,1)*params.alambda0_local;         
    Niters_batch = zeros(Ncfg,0);   
    iiter_total = 1;
    params.perform_final_iterations = false;
    flip_z_net = false;
    disp('reached')
    %% Go to get_pupil_matrix: they added two arguments for a function that takes one! This most likely causes the 'too many input arguments' error
    %% this is because you have the 3d OTF standalone, which has duplicate files of get_pupil_matrix!
    %% use the one with two arguments: it fits initialvalues_phasor
    %% Important: I replaced the get_otf file from vectorfitter with the file from the 3D_OTF_standalone file, as it is only viable for 2d-OTFs
    % disp(parameters)

    theta_init_batch = theta_batch;
    %disp('reached')
    Xpatch= loaded_input.Xpatch;
    Ypatch= loaded_input.Ypatch;
    Zpatch= loaded_input.Zpatch;
    params.Xpatch =  Xpatch;
    params.Ypatch =  Ypatch;
    params.Zpatch =  Zpatch;
    [theta_local,thetastore_local,localizations,localizations_with_outliers,meritstore_local,alambda_local,mu_allspots,dmudtheta_allspots,Niters_local,outliers] = ...
    local_update_otfmode(theta_init_batch,thetastore_local,theta_global,meritstore_local,alambda_local,spots_batch,roixy_batch,iiter_total,Niters_batch,flip_z_net,framelist(start+1:stop),ID(start+1:stop),OTF_3d_abs,params);
end                         
disp('success!')
% disp(mu)
% figure;
% imagesc(mu*(-1))
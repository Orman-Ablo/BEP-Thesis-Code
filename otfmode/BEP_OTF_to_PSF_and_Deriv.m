%% Step 1: Access proper files and variables
%Define input and output file locations
input_path ="C:\Users\Naam\Documents\BEP\Code\vectorfit-master\vectorfit-master\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_file.mat";
Example_file_spots = "C:\Users\Naam\Documents\BEP\Code\vectorfit-master\vectorfit-master\BEP Code Roman\Extra input data\transfer_3055388_files_4e4cf1e2\localization_segmentation_roi17_thr15p5_Data0001.mat";
output_folder = "C:\Users\Naam\Documents\BEP\Code\vectorfit-master\vectorfit-master\BEP Code Roman\Data_BEP_Roman\All_files_data\";
all_input_files = dir(input_path);
nr_of_datasets = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_datasets);

spots_information = dir(Example_file_spots);

%% set up parallel pool
nr_of_parallel_workers = 8;
if isempty(gcp('nocreate'))
    parpool('Threads', nr_of_parallel_workers);
end

%% Parameters: Use one of the parameters to initialize some values (e.g. file 1)
%% Input (Averaged PSFs and OTFs)


%% Load Parameters: OTF/PSF
dataset_name = all_input_files.name; 
foldername = all_input_files(1).folder;
input_filename = fullfile(foldername, dataset_name);
loaded_input = load(input_filename);
PSFmodelmatrix= loaded_input.PSFmodelmatrix;
Xpatch= loaded_input.Xpatch;
Ypatch= loaded_input.Ypatch;
Zpatch= loaded_input.Zpatch;
chisquarematrix= loaded_input.chisquarematrix;
numlocmatrix= loaded_input.numlocmatrix;
PSFdatamatrix= loaded_input.PSFdatamatrix;
OTF_3d = loaded_input.OTF_3d;
OTF_3d_abs = loaded_input.OTF_3d_abs;

%% Load Parameters: Spot Information
Spots_dataset = spots_information(1).name;
Spots_folder = spots_information(1).folder;
Spots_filename = fullfile(Spots_folder, Spots_dataset);
load_spots = load(Spots_filename);
localizations = load_spots.localizations;
theta_update = load_spots.theta_update;
outliers = load_spots.outliers;
allspots = load_spots.allspots;
roixy = load_spots.roixy;
mu = load_spots.mu;
params = load_spots.params;
Npatch = 4;%params.Npatch;



% define the limits of your PSF(spot), and based on this, return the proper
% OTF
 %% xy-limits
x_loc = localizations(1,1);% keep the row constant, but use a variable for the columns
y_loc = localizations(1,2);
z_loc = localizations(1,3);


for ix = 1:Npatch
    if (x_loc>Xpatch(ix))&&(x_loc<Xpatch(ix+1))
        x_patch = ix;
        break
    end
end

for iy = 1:Npatch
    if (y_loc>=Ypatch(iy))&&(y_loc<=Ypatch(iy+1))
        y_patch = iy;
        break
    end
end
disp(y_patch)
disp(x_patch)

% for zi = 1:params.mPSFz-1
%     if(zrange(zi) <= z_emit && z_emit<= zrange(zi+1))
%         z_patch = zi;
%         break
%     end
% end
OTF = OTF_3d_abs(:,:,:,x_patch,y_patch);
figure;
imagesc(OTF(:,:,10))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make sure all these parameters are present!
xemit = 0;%params.xemit;
yemit = 0;%params.yemit;
zemit = 0;%params.zemit;
params.xemit = xemit;
params.yemit = yemit;
params.zemit = zemit;
params.NA = 1.35;
params.Notf = params.Notfx;
params.Mroix = params.Mx;
params.Mroiy = params.My;
params.roisamplingdistance = params.pixelsize; % sampling distance in ROI
params.xroirange = params.roisamplingdistance*params.Mroix/2; % 1/2 of size PSF space in x
params.yroirange = params.roisamplingdistance*params.Mroiy/2; % 1/2 of size PSF space in y
params.fitmodel = 'xyz';
params.zstage = zemit;%% ask whther this makes sense
compders =1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PSF_Final,PSFderiv]= get_psfs_derivatives_otfmode(OTF,params,compders);
disp(PSF_Final) % why so many points are negative?
figure;
imagesc(PSF_Final)


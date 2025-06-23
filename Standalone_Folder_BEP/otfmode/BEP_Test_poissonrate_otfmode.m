%% Testing the parameters needed for, and effectiveness of, the poissonrate using 'get_psfs_derivatives_otfmode' to acquire the PSF and its derivatives



%input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_tukey_30.mat";
%input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_tukey_40.mat";
%input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_tukey_50.mat";
%input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_30.mat";
input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_40.mat";
%input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_50.mat";


all_input_files = dir(input_path);
nr_of_datasets = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_datasets);

Example_file_spots = "C:\Users\Naam\Documents\BEP\BEP Code Roman\input data\transfer_3055388_files_4e4cf1e2\localization_segmentation_roi17_thr15p5_Data0002.mat";

%% set up parallel pool
nr_of_parallel_workers = 8;
if isempty(gcp('nocreate'))
    parpool('Threads', nr_of_parallel_workers);
end




%% Parameters: Use one of the parameters to initialize some values (e.g. file 1)
%% Input (Averaged PSFs and OTFs)
dataset_name = all_input_files.name; 
foldername = all_input_files(1).folder;
input_filename = fullfile(foldername, dataset_name);
loaded_input = load(input_filename);
PSFmodelmatrix= loaded_input.PSFmodelmatrix;
PSFdatamatrix = loaded_input.PSFdatamatrix;
Xpatch= loaded_input.Xpatch;
Ypatch= loaded_input.Ypatch;
Zpatch= loaded_input.Zpatch;
chisquarematrix= loaded_input.chisquarematrix;
numlocmatrix= loaded_input.numlocmatrix;
OTF_3d = loaded_input.OTF_3d;
OTF_3d_abs = loaded_input.OTF_3d_abs;
Nph_all = loaded_input.Nph_all;
intensity = loaded_input.intensity;


all_spots_info = dir(Example_file_spots);
Spots_dataset = all_spots_info(1).name;
Spots_folder = all_spots_info(1).folder;
Spots_filename = fullfile(Spots_folder, Spots_dataset);
load_spots = load(Spots_filename);
localizations = load_spots.localizations;
theta_update = load_spots.theta_update;
outliers = load_spots.outliers;
allspots = load_spots.allspots;
roixy = load_spots.roixy;
mu = load_spots.mu;
params = load_spots.params;
% dipshow(allspots)
% diptruesize(2000)
% dipshow(mu)
% diptruesize(2000)

framelist = localizations(:,5);
ID = localizations(:,1);
Npatch = size(OTF_3d,5);%params.Npatch;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose the points to observe
jcfg = 250;%randi([1 1000],1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


theta_fit = load_spots.theta_update;

locindex_file = 1:size(allspots,4);
masklocs = ~ismember(locindex_file,outliers);
%% Makes Sure all the used variables are properly Seperated!
allspots = allspots(:,:,:,masklocs);
mu = mu(:,:,:,masklocs);
roixy = roixy(:,masklocs);
theta_update = theta_update(:,masklocs);



z_fit = [min(Zpatch) max(Zpatch)];
sep = find((localizations(:,4)>= z_fit(1)) & (localizations(:,4)<= z_fit(2)));
ID = ID(sep);
framelist = framelist(sep);
allspots = allspots(:,:,:,sep);
mu = mu(:,:,:,sep);

roixy = roixy(:,sep);
theta_update = theta_update(:,sep);
localizations = localizations(sep,:);





%% xy-limits
roix = roixy(1,jcfg);
roiy = roixy(2,jcfg);
x_loc = theta_update(1,jcfg)+roix*params.pixelsize;% keep the row constant, but use a variable for the columns
y_loc = theta_update(2,jcfg)+roiy*params.pixelsize;
% Localizations (xy) = thetaUpdate (xy)+roixy*pixelsize






x_patch = NaN;
y_patch = NaN;
for ix = 1:Npatch
    if (x_loc>=Xpatch(ix))&&(x_loc<=Xpatch(ix+1))
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
if isnan(x_patch) || isnan(y_patch)
    error('Patch index not found! x_loc=%.2f, y_loc=%.2f', x_loc, y_loc);
end
Nph_array = Nph_all(:,x_patch,y_patch);
OTF = OTF_3d(:,:,:,x_patch,y_patch);


% make sure all these parameters are present!
xemit = theta_update(1,jcfg);%params.xemit;
yemit = theta_update(2,jcfg);%params.yemit;
zemit = theta_update(3,jcfg);%params.zemit;
params.xemit = xemit;
params.yemit = yemit;
params.zemit = zemit;

% size(PSFmodelmatrix);
mPSFx = 85;
mPSFy = 85;
mPSFz = 53;
params.Mpsfx = mPSFx;
params.Mpsfy = mPSFy;
params.Mpsfz = mPSFz;
params.Mroix = 17;
params.Mroiy = 17;
params.Mx = params.Mroix;
params.My = params.Mroiy;
nOTF = size(OTF,1);
nOTFz = size(OTF,3);
params.Notfz = nOTFz;
params.fitmodel = 'xyz';
% params.psfsamplingdistance = 0.25*(params.lambda/params.NA/4); % sampling distance in PSF space
% params.psfsamplingdistancez = 0.1*(params.lambda/params.NA/4); % sampling distance in PSF space
 params.roisamplingdistance = params.pixelsize; % sampling distance in ROI
% params.psfxrange = params.psfsamplingdistance*params.Mpsfx/2; % 1/2 of size PSF space in x
% params.psfyrange = params.psfsamplingdistance*params.Mpsfy/2; % 1/2 of size PSF space in y
% params.psfzrange = params.psfsamplingdistancez*params.Mpsfz/2;
 params.Notf = nOTF; % #sampling points in OTF space, issue #1 form standard params
% params.spat_freq_cutoff = params.NA*2/params.lambda; % Lateral cutoff frequency of parameters
% n_imm = params.refimmnom;
params.xroirange = params.roisamplingdistance*params.Mroix/2; % 1/2 of size PSF space in x
params.yroirange = params.roisamplingdistance*params.Mroiy/2; % 1/2 of size PSF space in y
params.ztype = 'medium';
compders = 1;
params.compders = compders;
params.zstage = 0;
thetatry = theta_update(:,jcfg); %first spot
Ncfg = 2e3;

 %precompute the prechirpz
        Notf =params.Notf;
        xrange = params.xroirange;
        yrange = params.yroirange;
        NA = params.NA;
        lambda = params.lambda;
        Mx = params.Mx;
        My = params.My;
        
        % OTF space size and sampling (in diffraction units)
        OTFSize = 2.0; % OTF cutoff in xy-direction
        DxyOTF = 2*OTFSize/Notf;
        XYOTF = -OTFSize+DxyOTF/2:DxyOTF:OTFSize;
        [YOTF,XOTF] = meshgrid(XYOTF,XYOTF);
        
        
        % ROI size (in diffraction units)
        ROISizex = xrange*NA/lambda;
        ROISizey = yrange*NA/lambda;
        
        allOTFSize = [OTFSize OTFSize];
        allROISize = [ROISizex ROISizey];
        allNotf = [Notf Notf];
        allM = [Mx My];
        
        
        % calculate auxiliary vectors for chirpz
        % [Ax,Bx,Dx] = prechirpz(OTFSize,ROISizex,Notf,Mx);
        % [Ay,By,Dy] = prechirpz(OTFSize,ROISizey,Notf,My);
        [allA,allB,allD] = prechirpzn(allOTFSize,allROISize,allNotf,allM);
        params.allA = allA; % store the prechirps vectors, since it's the same computation for all uses in the script.
        params.allB = allB;
        params.allD = allD;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alambda_local =ones(Ncfg,1)*params.alambda0_local;
alambda = alambda_local(jcfg);
spot = allspots(:,:,1,jcfg);
mu_1 = mu(:,:,1,jcfg);


varfit = params.varfit;
theta_old = thetatry;

params.Zpatch =Zpatch;
params.pos = linspace(-params.Mpsfz/2,params.Mpsfz/2,params.Mpsfz);
intensity = mean(intensity, [2, 3]);

[mu,dmudtheta] = poissonrate_otfmode2(params,thetatry,OTF,intensity);


% figure;
% imagesc(dmudtheta(:,:,3)) %  PSFder scaled by Nph,
% 
% dipshow(dmudtheta(:,:,3))
% diptruesize(2000)
% colormap parula
% 
% figure;
% imagesc(dmudtheta(:,:,2))
% 
% dipshow(dmudtheta(:,:,2))
% diptruesize(2000)
% colormap parula
% 
% figure;
% imagesc(dmudtheta(:,:,1))
% 
% dipshow(dmudtheta(:,:,1))
% diptruesize(2000)
% colormap parula

[merit,grad,Hessian] = likelihood(params,spot,mu,dmudtheta,varfit);

[thetamin,thetamax] = thetalimits(params,thetatry);
thetaretry = (thetamax+thetamin)/2;

[theta_new,~] = thetaupdate(thetatry,thetamax,thetamin,thetaretry,grad,Hessian,alambda,params);


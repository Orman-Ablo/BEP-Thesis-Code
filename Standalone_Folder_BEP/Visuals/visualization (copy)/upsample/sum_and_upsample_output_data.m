close all
clear all

%path_output_data = '';

load(path_output_data);
%%
pixelsize = params.pixelsize;
imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;

% convert positions from nm to pixel units
theta_local = theta_full.local;
%theta_local = theta_update;
allpositions = [theta_local(2,:); theta_local(1,:)];
allpositions = allpositions / pixelsize;
z_positions = theta_local(3,:);

% Remove outliers
Ncfg_total = size(roixy,2);
no_outliers = setdiff(1:Ncfg_total,outliers_full{1});
%no_outliers = setdiff(1:Ncfg_total,outliers);
theta_local = theta_local(:,no_outliers);
allpositions = allpositions(:,no_outliers);
z_positions = z_positions(:,no_outliers);
mu = mu_full(:,:,:,no_outliers);
%mu = mu(:,:,:,no_outliers);
roixy = roixy(:,no_outliers);
allspots = allspots(:,:,:,no_outliers);

%% show ROI locations in FOV
figure
scatter3(roixy(1,:),roixy(2,:),theta_local(3,:));

%%
% select indices in specified region

upsfac = 10;
imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;
Npatch = 4;
%range = 120:220;
range = 1:170;
%zrange = [-300,-100];
zrange = [-400,-100];
%zrange = [200,1000];

figure
for xi = 1:Npatch
    for yi = 1:Npatch

        xrange = [floor((imgSizeX/Npatch)*(xi-1))+1,floor((imgSizeX/Npatch))*xi];
        yrange = [floor((imgSizeY/Npatch)*(yi-1))+1,floor((imgSizeY/Npatch))*yi];
        indices = find((roixy(1,:) >= xrange(1) & roixy(1,:) <= xrange(2) & roixy(2,:) >= yrange(1) & roixy(2,:) <= yrange(2) ...
            & z_positions > zrange(1) & z_positions < zrange(2) ));
        %indices = indices(indices<1e5);
        Ncfg = numel(indices);

        % select spots with these indices
        allspots_selected = allspots(:,:,:,indices);
        allpositions_selected = allpositions(:,indices);

        % feed spot stack to function for shift upsample and sum
        PSFsum_measured = sum_shift_ups_PSFs(allspots_selected,allpositions_selected,upsfac);

        % sum modeled PSFs
        
        % if model_psf
        %     patch_x = (xrange(2) + xrange(1))/2;
        %     patch_y = (yrange(2) + yrange(1))/2;
        %     fov_coordinates = get_fov_coordinates([0 0]',patch_x,patch_y,params_model);
        %     zernikeCoefficients = get_zernike_coefficients(fov_coordinates(1), fov_coordinates(2), theta_global_model, params_model);
        %     [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params_model,zernikeCoefficients);
        %     [PSFsum_modeled,~] = poissonrate(params_model,[0 0 (zrange(2) + zrange(1))/2 1e5 0],PupilMatrix,allzernikes,wavevector,wavevectorzimm);    
        % else
            allspots_modeled = mu;
            %allspots_modeled = mu_adapted_all;
            allspots_modeled = allspots_modeled(:,:,:,indices);
            PSFsum_modeled = sum_shift_ups_PSFs(allspots_modeled,allpositions_selected,upsfac);
        %end

        % inspect the final result
        subplot(Npatch,Npatch*2,(xi-1)*2*Npatch+2*yi-1);
        (xi-1)*2*Npatch+2*yi-1
        %scrsz = [1 1 1366 768];
        %set(gcf,'Position',round([0.15*scrsz(3) 0.50*scrsz(4) 0.8*scrsz(3) 0.35*scrsz(4)]));
        %imagesc(log(PSFsum_measured))
        %imagesc(PSFsum_measured)
        %imagesc(log(PSFsum_measured(range,range)))
        imagesc(PSFsum_measured(range,range))
        axis square
        axis off
        title(sprintf('(x=%i, y=%i) measured. #spots = %i',xi,yi,numel(indices)))
        ax = gca;
        ax.FontSize = 8;

        subplot(Npatch,Npatch*2,(xi-1)*2*Npatch+2*yi);
        (xi-1)*2*Npatch+2*yi
        %imagesc(log(PSFsum_modeled))
        %imagesc(PSFsum_modeled)
        %imagesc(log(PSFsum_modeled(range,range)))
        imagesc(PSFsum_modeled(range,range))
        title(sprintf('(x=%i, y=%i) modeled',xi,yi))
        axis square
        axis off
        ax = gca;
        ax.FontSize = 8;

    end
end


%% One PSF


path_beadaber = '';
path_smlmaber = '';

data_beadaber = load(path_beadaber);
data_smlmaber = load(path_smlmaber);
%%

params = data_smlmaber.params;
roixy = data_smlmaber.roixy;
allspots = data_smlmaber.allspots;
mu_beadaber = data_beadaber.mu;
mu_smlmaber = data_smlmaber.mu;
Ncfg = size(roixy,2);
outliers_beadaber = data_beadaber.outliers;
outliers_smlmaber = data_smlmaber.outliers;
outliers_combined = unique(cat(1,outliers_beadaber,outliers_smlmaber));
no_outliers = setdiff(1:Ncfg,outliers_combined);

roixy = roixy(:,no_outliers);
allspots = allspots(:,:,:,no_outliers);
mu_beadaber = mu_beadaber(:,:,:,no_outliers);
mu_smlmaber = mu_smlmaber(:,:,:,no_outliers);

imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;

theta_local_smlm = data_smlmaber.theta_update;
theta_local_smlm = theta_local_smlm(:,no_outliers);
%z_positions_smlm = theta_local_smlm(3,:);
allpositions_smlm = [theta_local_smlm(2,:); theta_local_smlm(1,:)];
allpositions_smlm = allpositions_smlm / params.pixelsize;

theta_local_bead = data_beadaber.theta_update;
theta_local_bead = theta_local_bead(:,no_outliers);
%z_positions_bead = theta_local_bead(3,:);
allpositions_bead = [theta_local_bead(2,:); theta_local_bead(1,:)];
allpositions_bead = allpositions_bead / params.pixelsize;

%%
yrange = [45,55]*params.imgSizeX/params.pixelsize; %pixels
xrange = [45,55]*params.imgSizeX/params.pixelsize; %pixels
%zrange = [-200, 200];
%zrange = [459,659]; %nm % 700 nm 
%zrange = [390,490]; %nm % 600 nm
%zrange = [-300,-200]; %A
%zrange = [-490,-390]; %B
%zrange = [480,680]; %nm % 400 nm 
%indices = find((roixy(1,:) >= xrange(1) & roixy(1,:) <= xrange(2) & roixy(2,:) >= yrange(1) & roixy(2,:) <= yrange(2)  & z_positions > zrange(1) & z_positions < zrange(2) ));
indices = find((roixy(1,:) >= xrange(1) & roixy(1,:) <= xrange(2) & roixy(2,:) >= yrange(1) & roixy(2,:) <= yrange(2) ));
%indices = indices(indices<1e5);
Ncfg = numel(indices)
%mean(z_positions(indices))
upsfac = 10;

% select spots with these indiced
allspots_selected = allspots(:,:,:,indices);
allpositions_selected_smlm = allpositions_smlm(:,indices);
allpositions_selected_bead = allpositions_bead(:,indices);

% feed spot stack to function for shift upsample and sum
PSFsum_measured_smlmpos = sum_shift_ups_PSFs(allspots_selected,allpositions_selected_smlm,upsfac);
PSFsum_measured_beadpos = sum_shift_ups_PSFs(allspots_selected,allpositions_selected_bead,upsfac);

% sum modeled PSFs
allspots_modeled_smlm = mu_smlmaber;
allspots_modeled_smlm = allspots_modeled_smlm(:,:,:,indices);
PSFsum_modeled_smlm = sum_shift_ups_PSFs(allspots_modeled_smlm,allpositions_selected_smlm,upsfac);

% allspots_modeled = mu_adapted_all;
% allspots_modeled = allspots_modeled(:,:,:,indices);
% PSFsum_modeledt2 = sum_shift_ups_PSFs(allspots_modeled,allpositions_selected,upsfac);
allspots_modeled_bead = mu_beadaber;
allspots_modeled_bead = allspots_modeled_bead(:,:,:,indices);
PSFsum_modeled_bead = sum_shift_ups_PSFs(allspots_modeled_bead,allpositions_selected_bead,upsfac);

% allspots_modeled = mu_adapted_all;
% allspots_modeled = allspots_modeled(:,:,:,indices);
% PSFsum_modeledt2 = sum_shift_ups_PSFs(allspots_modeled,allpositions_selected,upsfac);
% 
% allspots_modeled = mu_adapted_all_beads;
% allspots_modeled = allspots_modeled(:,:,:,indices);
% PSFsum_modeledbead = sum_shift_ups_PSFs(allspots_modeled,allpositions_selected,upsfac);


%imagesc(cat(2,PSFsum_measured,PSFsum_modeled,PSFsum_modeledt2,PSFsum_modeledbead))
imagesc(cat(2,PSFsum_measured,PSFsum_modeled,PSFsum_modeledbead))

% imagesc(cat(1,PSFsum_measured,PSFsum_modeled,PSFsum_modeledt2,PSFsum_modeledbead))
imagesc(cat(2,PSFsum_measured_smlmpos,PSFsum_modeled_smlm,PSFsum_measured_beadpos,PSFsum_modeled_bead))

axis equal
axis off
ax = gca;
ax.FontSize = 8;


% figure
% subplot(1,2,1)
% imagesc(PSFsum_measured)
% axis square
% axis off
% subplot(1,2,2)
% imagesc(PSFsum_modeled)
% axis square
% axis off
% ax = gca;
% ax.FontSize = 8;
%% Show averaged spots in 3D
figure
surf(PSFsum_measured);

%% Show all spots
measured_fitted = cat(2,allspots,mu);
dipshow(measured_fitted,'lin','name','simulated vs. fitted spots (without outliers)');
diptruesize(40000/params.Mx)
colormap parula

%%
PSFsum_normalized = 255*PSFsum_measured/max(PSFsum_measured,[],'all');
imwrite(uint8(PSFsum_normalized),'sumPSF.tiff')

%% model PSFs with arbitrary aberrations
Ncfg = size(allpositions,2);
if isempty(gcp('nocreate'))
    number_of_parallel_workers = feature('numcores');
    myCluster = parcluster('local');
    parpool(myCluster, number_of_parallel_workers);
end

theta_global = squeeze(theta_full.global(:,:,1));
%theta_global = squeeze(theta_global);
theta_global(12) = -18;
% theta_global = [       0
%          0
%          0
%          0
%    -0.4352
%    59.2322
%     0.7173
%     0.7955
%    -0.1002
%     4.7825
%    -2.9614
%   -17.9881
%    17.1018]; % Fitted aberrations by inversion from beads
% theta_global = theta_global.*[-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1]';
fov_coordinates = get_fov_coordinates(roixy,theta_local(1,:),theta_local(2,:),params);
fov_x = fov_coordinates(1,:);
fov_y = fov_coordinates(2,:);
zernikeCoefficients = get_zernike_coefficients(fov_x, fov_y, theta_global, params);

%mu_adapted_all_beads = zeros(size(allspots));
mu_adapted_all = zeros(size(allspots));
params.fitmode = 'PSF';
params.fitmodel = 'xyz';

parfor jcfg=1:Ncfg
    jcfg
    [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params,zernikeCoefficients(jcfg,:));
    [mu_adapted,~] = poissonrate(params,theta_local(:,jcfg),PupilMatrix,allzernikes,wavevector,wavevectorzimm);    
    %mu_adapted_all_beads(:,:,:,jcfg) = mu_adapted;
    mu_adapted_all(:,:,:,jcfg) = mu_adapted;
end


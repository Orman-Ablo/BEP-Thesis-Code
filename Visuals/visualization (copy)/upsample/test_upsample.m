% This script is for testing the function to sum shifted upsampled spots
% into an overall PSF.
%
% copyright Sjoerd Stallinga, TU Delft, 2023

close all
clear all

params = set_parameters_highresim_nanorulers;

Mx = params.Mx;
My = params.My;
Mz = params.Mz;

Ncfg = 1000;
allpositions = 5*(rand(2,1000)-0.5);
%allpositions = zeros(2,1000);

xy = allpositions*params.pixelsize;
bg = 10; % baclground level
N = 3000; % signal photon count


% Generate PSFs 
allspots = zeros(Mx,My,Mz,Ncfg);
for jcfg=1:Ncfg
    jcfg
    params.xemit = xy(1,jcfg);
    params.yemit = xy(2,jcfg);
    [wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = ...
    get_pupil_matrix(params);
    [FieldMatrix,FieldMatrixDerivatives] = ...
    get_field_matrix_derivatives(params,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
    [PSF,~] = get_psfs_derivatives(params,PupilMatrix,FieldMatrix,FieldMatrixDerivatives);
    PSF = N*PSF+bg;
    allspots(:,:,:,jcfg) = PSF;
end

% sig = 1.5; % spot width (pixel units)
% 
% % set up coordinate grid of ROI
% xx = (1:Mx)-(Mx+1)/2;
% yy = (1:My)-(My+1)/2;
% [Xg,Yg] = meshgrid(yy,xx);
% 
% allspots = zeros(Mx,My,Mz,Ncfg);
% for js = 1:Ncfg
%     js
%   centerpos = allpositions(:,js);
%   spot = exp(-((Xg-centerpos(1)).^2+(Yg-centerpos(2)).^2)/2/sig^2)/(2*pi*sig^2); % Gaussian
%   spot = bg+N*spot;
%   %allspots(:,:,js) = spot; % without noise
%   allspots(:,:,js) = poissrnd(spot); % with noise
% end

%%
dipshow(allspots)
diptruesize(4000)
colormap parula

%%
allpositions_new = [allpositions(2,:); allpositions(1,:)];

% feed spot stack to function for shift upsample and sum
upsfac = 10; % upsampling factor (must be integer)
PSFsum = sum_shift_ups_PSFs(allspots,allpositions_new,upsfac);

% inspect the final result
figure
imagesc(PSFsum)
axis square
axis off
colorbar
title('result sum/shift/upsample')

% check on signal strength
Nphotons_true = sum(allspots(:));
Nphotons_sum = sum(PSFsum(:));


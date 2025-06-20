% This script is for testing the function to sum shifted upsampled spots
% into an overall PSF.
%
% copyright Sjoerd Stallinga, TU Delft, 2023

close all
clear all

% input parameters
Nx = 9; % size ROIx
Ny = 9; % size ROIy
upsfac = 5; % upsampling factor (must be integer)
allpositions = [1.4,0.7;-1.1,0.2;1.4,-0.5]'; % true shift of spot center (in pixel units)
% allpositions = [2.4,1.7;2.4,1.7;2.4,1.7]'; % true shift of spot center (in pixel units)
% allpositions = [2.4,0;2.4,0;2.4,0]'; % true shift of spot center (in pixel units)
% allpositions = [0,1.7;0,1.7;0,1.7]'; % true shift of spot center (in pixel units)
%allpositions = [0,0;0,0;0,0]'; % true shift of spot center (in pixel units)
Nstack = size(allpositions,2); % number of spots
bg = 10; % baclground level
N = 3000; % signal photon count
sig = 1.5; % spot width (pixel units)

% set up coordinate grid of ROI
xx = (1:Nx)-(Nx+1)/2;
yy = (1:Ny)-(Ny+1)/2;
[Xg,Yg] = meshgrid(yy,xx);

% generate model spots
allspots = zeros(Nx,Ny,Nstack);
for js = 1:Nstack
  centerpos = allpositions(:,js);
  spot = exp(-((Xg-centerpos(1)).^2+(Yg-centerpos(2)).^2)/2/sig^2)/(2*pi*sig^2); % Gaussian
  spot = bg+N*spot;
%   allspots(:,:,js) = spot; % without noise
  allspots(:,:,js) = poissrnd(spot); % with noise
end

% feed spot stack to function for shift upsample and sum
PSFsum = sum_shift_ups_PSFs(allspots,allpositions,upsfac);

% inspect the final result
figure
scrsz = [1 1 1366 768];
set(gcf,'Position',round([0.15*scrsz(3) 0.50*scrsz(4) 0.8*scrsz(3) 0.35*scrsz(4)]));
subplot(1,3,1)
imagesc(allspots(:,:,1))
axis square
axis off
colorbar
title('1st spot')
subplot(1,3,2)
imagesc(allspots(:,:,Nstack))
axis square
axis off
colorbar
title('last spot')
subplot(1,3,3)
imagesc(PSFsum)
axis square
axis off
colorbar
title('result sum/shift/upsample')

% check on signal strength
Nphotons_true = sum(allspots(:));
Nphotons_sum = sum(PSFsum(:));


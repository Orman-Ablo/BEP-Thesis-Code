% This script generates synthetic single molecule data, adds noise to it
% and plots/saves the data.
%
% Autor: Isabel Droste, TU Delft, 2022

close all
clear all

%% Set parameters
%parameters = set_parameters_xy;
%parameters = set_parameters_xyz;
%parameters = set_parameters_xyz_aberrations;
parameters = set_parameters_xy_aberrations;
%parameters = set_parameters_zstack_bead_0215;

% Settings
Ncfg = 1;
save_data = false;
dest_file = 'N:\tnw\IST\QI\users\idroste\Data\synthetic_data\small_beads\bead_data.mat';

parameters.Ncfg = Ncfg;

% Overall rms aberration level
Wrms = 0.01; % lambda
%Wrms = 0;

% Total photon count and background range. A random number will be selected
% in this range. To select a fixed number, choose a range with the same
% start and end point.
Nrange = [2500, 2500];
Brange = [3,3];
%Nrange = [1,1]; % Total photon count
%Brange = [0,0]; % Number of background photons per pixel

% Include dipimage
addpath('C:\Users\idroste\Desktop\TUDelft\Code\diplib\share\DIPimage');
setenv('PATH',['C:\Users\idroste\Desktop\TUDelft\Code\diplib\bin',';',getenv('PATH')]);


%% Generate data

% Generate PSFs
% for random ground truth values of x, y, z, aberration coefficients, 
% depending on the fit mode.

fprintf('Generating ground truth instances\n')
[groundtruth,allspots] = generate_synthetic_beads(Ncfg,Wrms,Nrange,Brange,parameters);

% Add poisson noise.
allspots = poissrnd(allspots);

% Calculate the sum of the generated beads
sum_of_intensities = squeeze(sum(allspots,[1 2]));
total_sum = sum(sum_of_intensities);

% Generate roixy randomly (=pixel location in full image).
imgsizeX = 2048;
imgsizeY = 2048;
roixy = zeros(Ncfg,2);
roixy(:,1) = randi(imgsizeX,Ncfg,1);
roixy(:,2) = randi(imgsizeY,Ncfg,1);

parameters.FOV = imgsizeX; 

% Save data.
if save_data
    save(dest_file,'allspots','roixy','parameters','groundtruth');
end

%% Show data
dipshow(allspots(:,:,:,1),'lin');
diptruesize(6000)
colormap parula


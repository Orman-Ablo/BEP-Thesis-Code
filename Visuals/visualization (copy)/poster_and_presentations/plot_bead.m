% This script plots a z-stack bead image using diplib

% Settings
params = set_parameters_zstack_bead_pieter_isabel;
addpath('C:\Users\idroste\Desktop\TUDelft\Code\diplib\share\DIPimage');
setenv('PATH',['C:\Users\idroste\Desktop\TUDelft\Code\diplib\bin',';',getenv('PATH')]);

% Load data
file = '\bead\bead_pieter';
datain = ['read' file '.tif'];

InfoImage =  imfinfo(datain);
try
    allspots = double(LoadTiff16bit(datain, [1 length(InfoImage)]));
catch
    allspots = zeros(InfoImage(1).Width,InfoImage(1).Height,length(InfoImage));
    for j=1:length(InfoImage)
        allspots (:,:,j) = double(imread(datain, j));
    end
    disp('Warning: Used slower imread instead of LoadTiff')
end
allspots = (allspots-params.offset)/params.gain;
allspots(allspots<=0) = 1e-3;

% Plot bead
dipshow(allspots,'lin');
diptruesize(800)
colormap hot


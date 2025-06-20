%% Render final image

%% Load data

%path_before_dc = '/home/idroste/Desktop/TUDelft/Data/Hertenlab/localization/fitted_aberrations/localizations_combined.mat';
%path_fitted_aberrations = '/home/idroste/Desktop/TUDelft/Data/Hertenlab/localization_drift_corrected/localizations_fitted_aberrations.mat';
%path_no_aberrations = '/home/idroste/Desktop/TUDelft/Data/Hertenlab/localization_drift_corrected/localizations_no_aberrations.mat';

path_fitted_aberrations = '/home/idroste/Desktop/TUDelft/Data/Lidke/data_analysis/DATA1/Data_wo_Astigmatism/localization/roi9/localization_combined.mat';
path_no_aberrations = '/home/idroste/Desktop/TUDelft/Data/Lidke/data_analysis/DATA1/Data_wo_Astigmatism/localization/roi9_noaber/localization_combined.mat';

%input_before_dc = load(path_before_dc);
input_fitted_aberrations = load(path_fitted_aberrations);
input_no_aberrations = load(path_no_aberrations);

%localizations_before_dc = input_before_dc.localizations;
localizations_fitted_aberrations = input_fitted_aberrations.localizations;
localizations_no_aberrations = input_no_aberrations.localizations;

params = input_fitted_aberrations.params;

%% Load data
clear all
close all

path_1 = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\Data_w_Astigmatism\localization_drift_corrected\localization_fitted_aberrations.mat';
path_2 = 'U:\NATaberrations\ExperimentalData\Lidke\raw_data\UNM_Data1\DATA1\Data_driftc_sml.mat';

input_1 = load(path_1);
input_2 = load(path_2);

localizations_1 = input_1.localizations;
x1 = localizations_1(:,2);
y1 = localizations_1(:,3);
z1 = localizations_1(:,4);

% open SMAP data
x2 = input_2.saveloc.loc.xnm;
y2 = input_2.saveloc.loc.ynm;
z2 = input_2.saveloc.loc.znm;

params = input_1.params;

%%
params.pixelsize = 95;
params.imgSizeX = 1024;
params.imgSizeY = 1024;
sizeX = double(params.pixelsize)*double(params.imgSizeX);
sizeY = double(params.pixelsize)*double(params.imgSizeY);

x1 = localizations_1(:,2);
y1 = localizations_1(:,3);
z1 = localizations_1(:,4);
x2 = localizations_2(:,2);
y2 = localizations_2(:,3);
z2 = localizations_2(:,4);

%rangex = [0 sizeX];
%rangey = [0 sizeY];
%rangex = [0.6*sizeX 0.8*sizeX];
%rangey = [0.2*sizeY 0.4*sizeY];
%rangex = [0.15*sizeX 0.25*sizeX];
%rangey = [0.25*sizeY 0.35*sizeY];
%rangex = [0.0*sizeX 0.15*sizeX];
%rangey = [0.0*sizeY 0.15*sizeY];
rangex = [0.0*sizeX 1.0*sizeX];
rangey = [0.0*sizeY 1.0*sizeY];

pixelsize = 5; % Should be ~half the size of the crlb. The smaller, the large the rendered image becomes.
alpha = 0.0;
beta = 0.0;
sigma = 2.0;
cmax = 3.0;
img_1 = renderprojection(x1,y1,z1,alpha,beta,rangex,rangey,pixelsize,sigma,false,cmax);
img_2 = renderprojection(x2,y2,z2,alpha,beta,rangex,rangey,pixelsize,sigma,false,cmax);

img_fa = renderprojection(x1,y1,z1,alpha,beta,rangex,rangey,pixelsize,sigma,false,cmax);
img_na = renderprojection(x2,y2,z2,alpha,beta,rangex,rangey,pixelsize,sigma,false,cmax);

%%
cmax = 1;
fig = figure(); 
ax(1) = axes('Units','normalized','Position', [ .0 .1 .7 .7]);
ax(2) = axes('Units','normalized','Position', [ .4 .1 .7 .7]);
%ax(3) = axes('Units','normalized','Position', [ .6 .3 .5 .5]);
imshow(img_fa,[0, cmax], 'Parent', ax(1))
imshow(img_na,[0, cmax], 'Parent', ax(2))
%imshow(img_na,[0, cmax], 'Parent', ax(3))
title(ax(1), 'Fitted aberrations')
title(ax(2), 'Constant aberrations')
%title(ax(3), 'No aberrations')
ax(1) = axes('Units','normalized','Position', [ .0 .3 .5 .5]);
ax(2) = axes('Units','normalized','Position', [ .3 .3 .5 .5]);
imshow(img_1,[0, cmax], 'Parent', ax(1))
imshow(0.7*img_2,[0, cmax], 'Parent', ax(2))
title(ax(1), 'Data 1')
title(ax(2), 'Data 2')
colormap hot
linkaxes(ax)
axis on

%%

[C,ia,ib] = intersect(localizations_fitted_aberrations(:,1),localizations_no_aberrations(:,1));

x_fa = localizations_fitted_aberrations(ia,2);
y_fa = localizations_fitted_aberrations(ia,3);
x_na = localizations_no_aberrations(ib,2);
y_na = localizations_no_aberrations(ib,3);

diff_x = x_fa - x_na;
diff_y = y_fa - y_na;

figure
hist(diff_y,10000)
xlim([-3 7])
xlabel('nm')
title('Histogram x-difference')
ax=gca;
ax.FontSize = 16;

figure
hist(diff_y,10000)
xlim([-3 7])
xlabel('nm')
title('Histogram y-difference')

%%

figure
plot(input_2.saveloc.file.driftinfo.xy.dyt)
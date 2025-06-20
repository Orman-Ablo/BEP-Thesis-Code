%% Render final image

%% Combine localizations into one array

% input_folder = '/home/idroste/Desktop/TUDelft/Data/yutong_nanorulers/data_analysis/localization/localization_*';
% 
% all_input_files = dir(input_folder);
% nr_of_files = size(all_input_files,1);
% fprintf('Found %i files\n',nr_of_files);
% 
% localizations = [];
% localizations_with_outliers = [];
% theta_init = [];
% theta = [];
% outliers = [];
% roixy = [];
% allspots = [];
% mu = [];
% Ncfg_total = 0;
% 
% for file_i=1:nr_of_files
% 
%     % Load data
%     filename = all_input_files(file_i).name; 
%     foldername = all_input_files(file_i).folder;
%     path_input_data_full = fullfile(foldername, filename);    
%     fprintf('Input data = %s\n',path_input_data_full);
%     loaded_data = load(path_input_data_full);
% 
%     localizations_temp = loaded_data.localizations;
%     localizations_with_outliers_temp = loaded_data.localizations_with_outliers;
%     theta_init_temp = loaded_data.theta_init;
%     outliers_temp = loaded_data.outliers;
%     roixy_temp = loaded_data.roixy;
%     allspots_temp = loaded_data.allspots;
%     mu_temp = loaded_data.mu;
%     theta_temp = loaded_data.theta_update;
%     params = loaded_data.params;
% 
%     localizations = cat(1,localizations,localizations_temp);
%     localizations_with_outliers = cat(1,localizations_with_outliers,localizations_with_outliers_temp);
%     theta_init = cat(2,theta_init,theta_init_temp);
%     outliers_temp_adjusted = outliers_temp + Ncfg_total;
%     outliers = cat(1,outliers,outliers_temp_adjusted);
%     Ncfg_total = Ncfg_total + size(roixy_temp,2);
%     roixy = cat(2,roixy,roixy_temp);
%     allspots = cat(4,allspots,allspots_temp);
%     mu = cat(4,allspots,allspots_temp);
%     theta = cat(2,theta,theta_temp);
% end
% 
% params = loaded_data.params;
% 
% % path_input_data = 'C:\Users\idroste\Desktop\TUDelft\Data\NAT_data\Fu_data_temp\drift_correction\RCC\drift_corrected_localizations_all.mat';
% % load(path_input_data)
% 
% %%
% output_folder = '/home/idroste/Desktop/TUDelft/Data/yutong_nanorulers/data_analysis/localization_combined/localization_combined_ROI7_NR80R_1114_Oil100x1-5_Col200_647P200-100_80ms_EM100_10000f_FL001.mat';
% save(output_folder,'theta','localizations','localizations_with_outliers','theta_init','outliers','roixy','allspots','mu','params','-v7')
% 
% %% Apply drift correction
% 
% %drift = importdata('/home/idroste/Desktop/TUDelft/Data/yutong_nanorulers/data_analysis/localization_drift_corrected/drift_output_yutong.txt');
% drift = importdata('/home/idroste/Desktop/TUDelft/Code/drift_correction/DME/drift-estimation-master/drift_output_highresim002.txt');
% 
% 
% localizations_before_dc = localizations;
% 
% %% 
% path = '/home/idroste/Desktop/TUDelft/Data/Jungmann/data_analysis/localization_drift_corrected/localization_dc_100fpb_1e7spots.mat';
% driftdata = load(path);
% drift = driftdata.drift;
% %%
% localizations_drift_corrected = localizations;
% 
% sz = size(localizations_drift_corrected);
% for ii=1:sz(1)
%     t = localizations_drift_corrected(ii,5);
%     localizations_drift_corrected(ii,2) = localizations_drift_corrected(ii,2) - drift(t,1);
%     localizations_drift_corrected(ii,3) = localizations_drift_corrected(ii,3) - drift(t,2);
%     %localizations_drift_corrected(ii,4) = localizations_drift_corrected(ii,4) - drift(t,3);
% end
% 
% %%
% localizations = localizations_drift_corrected;
% %%
% localizations_drift_corrected = localizations;
% 
% sz = size(localizations_drift_corrected);
% for ii=1:sz(1)
%     t = localizations_drift_corrected(ii,5);
%     localizations_drift_corrected(ii,2) = localizations_drift_corrected(ii,2) - drift(t,1);
%     localizations_drift_corrected(ii,3) = localizations_drift_corrected(ii,3) - drift(t,2);
% end
% %%
% localizations = localizations_drift_corrected;
% 
% %% Apply drift correction
% path = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Jungmann\data_analysis\find_localization_bug\no_aberrations\no_OTF\gpu_lowaccuracy\localization_drift_corrected\localization_dc_100fpb_1e7spots_drift_used.mat';
% driftdata = load(path);
% drift = driftdata.drift;
% 
% %%
% localizations_drift_corrected = localizations;
% 
% sz = size(localizations_drift_corrected);
% for ii=1:sz(1)
%     t = localizations_drift_corrected(ii,5);
%     localizations_drift_corrected(ii,2) = localizations_drift_corrected(ii,2) - drift(t,1);
%     localizations_drift_corrected(ii,3) = localizations_drift_corrected(ii,3) - drift(t,2);
% end
% 
% 
% %%
% localizations = localizations_drift_corrected;
% 
% %%
% 
% idx_highNph = find(localizations(:,9)>2500);
% %idx_highNph = find(localizations(:,9)>2500 & localizations(:,6)<5);
% %idx_highNph = find(localizations(:,6)<5);
% localizations_highNph = localizations(idx_highNph,:);
% 
% %%
% figure
% histogram(localizations(:,9))
% xlim([-1e3 1e4])
% 
% figure
% histogram(localizations_highNph(:,9))
% xlim([-1e3 1e4])


%%
%params = set_parameters_Jungmann;
%params = set_parameters_Lidke_2D;
params = set_parameters_nanorulers_hri_2D;


%x = coords_connect(:,1);
%y = coords_connect(:,2);
%z = coords_connect(:,3);

% x1 = localizations.y*params.pixelsize;
% y1 = localizations.x*params.pixelsize;
% z1 = zeros(size(x));
Nmax = size(localizations,1);%5e5;
% x2 = localizations(1:Nmax,2);
% y2 = localizations(1:Nmax,3);
% z2 = localizations(1:Nmax,4);

x1 = localizations(1:Nmax,2);
y1 = localizations(1:Nmax,3);
z1 = localizations(1:Nmax,4);


% x2 = localizations_before_dc(1:Nmax,2);
% y2 = localizations_before_dc(1:Nmax,3);
% z2 = localizations_before_dc(1:Nmax,4);

%Nmin = 1;%5.8281e+07:
%Nmax = 5.8281e+07;%116561505;

%x = localizations_drift_corrected(1:Nmax,2);
%y = localizations_drift_corrected(1:Nmax,3);
%z = localizations_drift_corrected(1:Nmax,4);

%x = localizations_drift_corrected_matlab(:,2);
%y = localizations_drift_corrected_matlab(:,3);
%z = localizations_drift_corrected_matlab(:,4);

%x = localizations_drift_corrected(:,2);
%y = localizations_drift_corrected(:,3);
%z = localizations_drift_corrected(:,4);

% x = saveloc.loc.xnm;
% y = saveloc.loc.ynm;
% z = zeros(size(x));

% x = localizations(:,2);
% y = localizations(:,3);
% z = localizations(:,4);

% y = localizations.x*params.pixelsize;
% x = localizations.y*params.pixelsize;
% z = zeros(size(x));
 
%y = localizations_drift_corrected.x*params.pixelsize + params.pixelsize;
%x = localizations_drift_corrected.y*params.pixelsize + params.pixelsize;
%z = zeros(size(x));

% x = localizations_vectorfit(:,2);
% y = localizations_vectorfit(:,3);
% z = localizations_vectorfit(:,4);

%x = localizations_with_outliers(:,2);
%y = localizations_with_outliers(:,3);
%z = localizations_with_outliers(:,4);

%x = localizations_highNph(:,2);
%y = localizations_highNph(:,3);
%z = localizations_highNph(:,4);

sizeX = params.pixelsize*params.imgSizeX;
sizeY = params.pixelsize*params.imgSizeY;

rangex = [0.0*sizeX 1.0*sizeX];
rangey = [0.0*sizeY 1.0*sizeY];

%rangex = [0.65*params.pixelsize*params.imgSizeY 0.69*params.pixelsize*params.imgSizeX];
%rangey = [0.48*params.pixelsize*params.imgSizeY 0.53*params.pixelsize*params.imgSizeY];
%rangey = [0 1000];
pixelsize = 5%mean(localizations(:,6:7),'all')/2 % Should be ~half the size of the crlb. The smaller, the large the rendered image becomes.
alpha = 0.0;
beta = 0.0;
sigma = 2.0;
cmax = 4.0;
rendered_img_1 = renderprojection(x1,y1,z1,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax); %shift 1 pixel for Lidke data
%rendered_img_2 = renderprojection(x2,y2,z2,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax);


%% Compare two images
% 
cmax = 4.0;
fig = figure(); 
subplot(1,2,1)
imshow(0.8*rendered_img_1,[0, cmax])
axis on
title('With aberrations')
subplot(1,2,2)
imshow(0.8*rendered_img_2,[0, cmax])
axis on
title('Without aberrations')
colormap hot
linkaxes([subplot(1,2,1),subplot(1,2,2)])

%cmax = 1.5;
%rendered_img_1 = renderprojection(x,y,z,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax);

%% Compare two images

cmax = 4.0;
fig = figure(); 
ax(1) = axes('Units','normalized','Position', [ .0 .3 .5 .6]);
ax(2) = axes('Units','normalized','Position', [ .5 .3 .5 .6]);
imshow(1*rendered_img_1,[0, cmax], 'Parent', ax(1))
imshow(1*rendered_img_2,[0, cmax], 'Parent', ax(2))
title(ax(1), 'Picasso')
title(ax(2), 'Vectorfit')
colormap hot
linkaxes(ax)
axis on

%%
%x = localizations_ourdata(:,2);
%y = localizations_ourdata(:,3);
%z = localizations_ourdata(:,4);

% params.pixelsize = 130;
% params.imgSizeX = 1024;
% params.imgSizeY = 1024;

%x = localizations_drift_corrected.x*params.pixelsize;
%y = localizations_drift_corrected.y*params.pixelsize;
%z = zeros(size(x));

%x = coords_connect(:,1);
%y = coords_connect(:,2);
%z = coords_connect(:,3);

%x = localizations_drift_corrected(1:Nmax,2);
%y = localizations_drift_corrected(1:Nmax,3);
%z = localizations_drift_corrected(1:Nmax,4);

%params.pixelsize = 130;
%params.imgSizeX = 1024;
%params.imgSizeY = 1024;

x = localizations_drift_corrected.x*params.pixelsize;
y = localizations_drift_corrected.y*params.pixelsize;
z = zeros(size(x));

% start = 1;
% stop = 13127995;%floor(size(x,1)/2);
% x = x(start:stop);
% y = y(start:stop);
% z = z(start:stop);

idx = find(localizations.frame<4074);
x = x(idx);
y = y(idx);
z = z(idx);

% x = xall;
% y = yall;
% z = zeros(size(x));

%x = localizations_no_aberrations_drift_corrected(:,2);
%y = localizations_no_aberrations_drift_corrected(:,3);
%z = localizations_no_aberrations_drift_corrected(:,4);


%x = localizations.x*params.pixelsize;
%y = localizations.y*params.pixelsize;
%z = zeros(size(x));

sizeX = double(params.pixelsize)*double(params.imgSizeX);
sizeY = double(params.pixelsize)*double(params.imgSizeY);

%rangex = [0 sizeX];
%rangey = [0 sizeY];
rangex = [0.4*params.pixelsize*params.imgSizeY 0.5*params.pixelsize*params.imgSizeX];
rangey = [0.4*params.pixelsize*params.imgSizeY 0.5*params.pixelsize*params.imgSizeY];
%rangey = [0 1000];
pixelsize = 1.7;%3.4 % Should be ~half the size of the crlb. The smaller, the large the rendered image becomes.
alpha = 0.0;
beta = 0.0;
sigma = 2.0;
cmax = 0.4;
image_rendered = renderprojection(y,x,z,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax);



%%
params.pixelsize = 95;
params.imgSizeX = 1024;
params.imgSizeY = 1024;

x = saveloc.loc.xnm;
y = saveloc.loc.ynm;
z = zeros(size(x));

sizeX = double(params.pixelsize)*double(params.imgSizeX);
sizeY = double(params.pixelsize)*double(params.imgSizeY);

rangex = [0 sizeX];
rangey = [0 sizeY];

pixelsize = 5;
alpha = 0.0;
beta = 0.0;
sigma = 2.0;
cmax = 2;
image_rendered = renderprojection(x,y,z,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax);

%% zoomed in images
rangex = [0.1*sizeX 0.35*sizeX];
rangey = [0.1*sizeY 0.35*sizeY];
pixelsize = 7;
image_rendered = renderprojection(x,y,z,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax);

%%
rangex = [0.15*sizeX 0.25*sizeX];
rangey = [0.25*sizeY 0.35*sizeY];
pixelsize = 5;
cmax = 4.0;
image_rendered = renderprojection(x,y,z,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax);

%%
rangex = [0.1*sizeX 0.3*sizeX];
rangey = [0.0*sizeY 0.2*sizeY];
pixelsize = 5;
cmax = 2.0;
image_rendered = renderprojection(x,y,z,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax);
%%
image_rendered = renderprojection(x_bc,y_bc,z,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax);

%%
rangex = [0.6*sizeX 0.8*sizeX];
rangey = [0.2*sizeY 0.4*sizeY];
pixelsize = 5;
sigma = 1.5;
cmax = 3.0;
image_rendered = renderprojection(x,y,z,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax);
%image_rendered = renderprojection(x_bc,y_bc,z,alpha,beta,rangex,rangey,pixelsize,sigma,true,cmax);

%%
figure
surf(image_rendered(1:20:13312,1:20:13312))
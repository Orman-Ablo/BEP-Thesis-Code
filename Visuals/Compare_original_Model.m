folder = "C:\Users\Naam\Documents\BEP\data\second_fit\";
folder_original = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Extra input data\transfer_3055388_files_4e4cf1e2\";
Patches = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_40.mat";
load(Patches,'Xpatch','Ypatch')


% Below is an older version that computed the combined locations, then
% poltted. The new version combines files preemptively.



% all_input_files = dir(folder);
% % Remove current and parent directories (show up as files across the
% % directory
% all_input_files = all_input_files(~[all_input_files.isdir]);
% 
% nr_of_files = size(all_input_files,1);
% fprintf('Found %i files\n',nr_of_files);
% 
% xstore1 = [];
% ystore1 =[];
% zstore1 = [];
% thetax1 = [];
% thetay1 =[];
% thetaz1 = [];
% for dataset_i = 1:100
%     dataset_name = all_input_files(dataset_i).name; 
%     foldername = all_input_files(dataset_i).folder;
%     input_filename = fullfile(foldername, dataset_name);
%     fprintf('Input file = %s\n',dataset_name);
%     loaded_input = load(input_filename);
%     localizations = loaded_input.localizations;
%     outliers = loaded_input.outliers;
% 
%     theta_update = loaded_input.theta_update;
%     filter = loaded_input.merit;
%     locindex_file = 1:size(filter,1);
%     masklocs = setdiff(locindex_file,outliers);
%     x_theta = theta_update(1,masklocs);
%     y_theta = theta_update(2,masklocs);
%     z_theta = theta_update(3,masklocs);
% 
%     thetax1 = [thetax1 x_theta];
%     thetay1 = [thetay1 y_theta];
%     thetaz1 = [thetaz1 z_theta];
%     params = loaded_input.params;
%     xtemp = localizations(:,2);
%     ytemp = localizations(:,3);
%     ztemp  = localizations(:,4);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % roixy = loaded_input.roixy;
%     % roixy_temp = roixy;
%     % locindex_file = 1:size(roixy_temp,2);
%     % masklocs = setdiff(locindex_file,outliers);
%     % roixy_add = roixy_temp(:,masklocs);
%     % roixy1 = [roixy1 roixy_add];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%     xstore1 = [xstore1; xtemp];
%     ystore1 = [ystore1; ytemp];
%     zstore1 = [zstore1; ztemp];
% end
% 
% 
% 
% 
% orig_files = dir(folder_original);
% all_input_files = orig_files(~[orig_files.isdir]);
% 
% nr_of_files = size(all_input_files,1);
% fprintf('Found %i files\n',nr_of_files);
% 
% xstore2 = [];
% ystore2 =[];
% zstore2 = [];
% roixy2 = [];
% thetax2 = [];
% thetay2 =[];
% thetaz2 = [];
% for dataset_i = 1:100
%     dataset_name = all_input_files(dataset_i).name; 
%     foldername = all_input_files(dataset_i).folder;
%     input_filename = fullfile(foldername, dataset_name);
%     fprintf('Input file = %s\n',dataset_name);
%     loaded_input = load(input_filename);
%     localizations = loaded_input.localizations;
%     outliers = loaded_input.outliers;
% 
%     filter = loaded_input.merit;
% 
% 
%     theta_update = loaded_input.theta_update;
% 
%     locindex_file = 1:size(filter,2);
%     masklocs = setdiff(locindex_file,outliers);
%     x_theta = theta_update(1,masklocs);
%     y_theta = theta_update(2,masklocs);
%     z_theta = theta_update(3,masklocs);
% 
%     thetax2 = [thetax2; x_theta;];
%     thetay2 = [thetay2; y_theta;];
%     thetaz2 = [thetaz2; z_theta;];
% 
% 
% 
%     params = loaded_input.params;
%     xtemp = localizations(:,2);
%     ytemp = localizations(:,3);
%     ztemp  = localizations(:,4);
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     roixy = loaded_input.roixy;
%     roixy_temp = roixy;
%     locindex_file = 1:size(roixy_temp,2);
%     masklocs = setdiff(locindex_file, outliers);
%     roixy_add = roixy_temp(:,masklocs);
%     roixy2 = [roixy2 roixy_add];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     xstore2 = [xstore2; xtemp];
%     ystore2 = [ystore2; ytemp];
%     zstore2 = [zstore2; ztemp];
% end
% % localizations = 
% locations1 = [xstore1 ystore1 zstore1];
% 
% 
% theta1 = [thetax1; thetay1; thetaz1;];
% theta2 = [thetax2; thetay2; thetaz2;];
% 
% disp(min(locations1(:,1)))
% disp(max(locations1(:,1)))
% locations2 = [xstore2 ystore2 zstore2];

input1 =  "C:\Users\Naam\Documents\BEP\BEP Code Roman\data\Combined\localization_combined_initial_fit.mat";
load_initial = load(input1);
input2 = "C:\Users\Naam\Documents\BEP\BEP Code Roman\data\Combined\localization_combined_second_fit_scaled.mat";
load_second = load(input2);

locations1 = load_initial.localizations_drift_corrected;
locations1 = locations1(:,2:4);
locations2 = load_second.localizations_drift_corrected;
locations2 = locations2(:,2:4);
theta1 = load_initial.theta;
theta2 = load_second.theta;

ncolors = 19;
colorrange = [-500,400];

xlim = [0.0 1.0];
ylim = [0.0 1.0];
pixelsize = 50;

sigma = 2.0;
clim = 0.95; % the brightest (1-clim) fraction of pixels gets clipped 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %choose image slice 
pos1 = (locations1(:,1:2)-theta1(1:2,:)')/params.pixelsize;
pos2 = (locations2(:,1:2)-theta2(1:2,:)')/params.pixelsize;

pix_coordsx = [200 700]; % from 0 to 1024 (imgsizex)
pix_coordsy = [200 700]; % from 0 to 1024 (imgsizey)
%locations = distance in 

indices1 = find( (pix_coordsx(1) <= pos1(:,1)) & (pos1(:,1) <= pix_coordsx(2)) & (pix_coordsy(1) <= pos1(:,2)) & (pos1(:,2) <= pix_coordsy(2)));
indices2 = find( (pix_coordsx(1) <= pos2(:,1)) & (pos2(:,1) <= pix_coordsx(2)) & (pix_coordsy(1) <= pos2(:,2)) & (pos2(:,2) <= pix_coordsy(2)));



locations1 = locations1(indices1,:);
locations2 = locations2(indices2,:);

disp(min(locations1));
disp(max(locations2));

xlim = pix_coordsx/1024; % normalize between 0 and 1
ylim = pix_coordsy/1024; % normalize between 0 and 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale_bar = 5000; 
 
scale_bar_length_px = round(scale_bar / pixelsize);


figure;
img1 =render_3D_func(locations1,ncolors,colorrange,xlim,ylim,pixelsize,sigma,clim,params);
% set position of bar 
y_position = size(img1, 1) - 20;
x_position = 0 + scale_bar_length_px +2 ; 

% scale bar
rectangle('Position', [x_position, y_position, scale_bar_length_px, 5], ...
          'FaceColor', 'w', 'EdgeColor', 'none');

% scale bar length
text(x_position + scale_bar_length_px/2, y_position - 20, ...
     '5 μm', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', 22);
title('Intitial Fit');
ax = gca;
ax.FontSize = 22;

figure;
img2 =render_3D_func(locations2,ncolors,colorrange,xlim,ylim,pixelsize,sigma,clim,params);
% set position of bar 
y_position = size(img2, 1) - 20;
x_position = 0 + scale_bar_length_px +2 ; 

% scale bar
rectangle('Position', [x_position, y_position, scale_bar_length_px, 5], ...
          'FaceColor', 'w', 'EdgeColor', 'none');

% scale bar length
text(x_position + scale_bar_length_px/2, y_position - 20, ...
     '5 μm', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', 22);
title('OTF Fit');
ax = gca;
ax.FontSize = 22;



gray1 = mat2gray(rgb2gray(img1));
gray2 = mat2gray(rgb2gray(img2));

%Create RGB image from Two datasets
combined_img = zeros(size(img1));
combined_img(:,:,1) = gray2;  % red channel
combined_img(:,:,2) = gray1;  % green channel
combined_img(:,:,3) = 0;      % blue channel

% Show
figure;
imshow(combined_img);



hold on
% set position of bar 
y_position = size(combined_img, 1) - 20;
x_position = 0 + scale_bar_length_px +2 ; 

% scale bar
rectangle('Position', [x_position, y_position, scale_bar_length_px, 5], ...
          'FaceColor', 'w', 'EdgeColor', 'none');

% scale bar length
text(x_position + scale_bar_length_px/2, y_position - 20, ...
     '5 μm', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', 22);
title('Comparison of OTF (red) vs First (Green) Fit');

ax = gca;
ax.FontSize = 22;


% for checking image limits
% check_file = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store.mat";
% load(check_file,'Xpatch','Ypatch','Zpatch')
% Npatch = 4;
% for ix = 1:Npatch
%     for iy = 1:Npatch
%          positionmask = (locations2(:,1)>=Xpatch(ix))&(locations2(:,1)<=Xpatch(ix+1))&...
%                        (locations2(:,2)>=Ypatch(iy))&(locations2(:,2)<=Ypatch(iy+1));
%          locations_z = locations2(positionmask,3);
%          max_z = max(locations_z);
%          min_z = min(locations_z);
%          sprintf('for patch (%i,%i), the minimum is %i, and the maximum is %i',ix,iy,min_z,max_z)
% 
%     end
% end

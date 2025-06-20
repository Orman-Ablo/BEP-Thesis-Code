% Make contour plot of aberration surfaces and compare to bead surfaces

clear all
close all

%%

% load('/home/idroste/Desktop/TUDelft/Data/Fu_fig4/data_analysis/fitted_aberrations/009_fitted_aberrations_segmentation_roi17_thr25_NUP96_SNP647_astigmatism_3D_1608_10ms_hama_mm_1800mW_3_MMStack_Default_009.ome.mat');
% load('/home/idroste/Desktop/TUDelft/Data/Fu_fig4/raw_data/data_paper/aber_map.mat');

%load('U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_fig4\data_analysis\fitted_aberrations\009_fitted_aberrations_segmentation_roi17_thr25_NUP96_SNP647_astigmatism_3D_1608_10ms_hama_mm_1800mW_3_MMStack_Default_009.ome.mat');
%load('U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_fig4\raw_data\data_paper\aber_map.mat');


%%
theta.global = theta_full.global;

fitmodel = params.fitmodel;
xn = linspace(-1,1,100);
yn = linspace(-1,1,100);

[Xn,Yn] = meshgrid(xn,yn);


nr_of_zernikes = size(params.aberrations,1);
orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernikes,1);

for jzer = 1:nr_of_zernikes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

% Zernike surfaces
if params.groundtruth_exists
    zernikeSurfaces_groundtruth = get_zernike_coefficients(Xn,Yn,groundtruth.global,params);
end
zernikeSurfaces_fitted = get_zernike_coefficients(Xn,Yn,theta.global,params);

% Convert aberrations from nm to mlambda
if params.groundtruth_exists
    zernikeSurfaces_groundtruth_normalized = 1e3*zernikeSurfaces_groundtruth/params.lambda;
end
zernikeSurfaces_fitted_normalized = 1e3*zernikeSurfaces_fitted/params.lambda;

if strcmp(fitmodel,'xy-gamma')
    if params.groundtruth_exists
        zmin = min(cat(1,zernikeSurfaces_groundtruth_normalized,zernikeSurfaces_fitted_normalized),[],'all');
        zmax = max(cat(1,zernikeSurfaces_groundtruth_normalized,zernikeSurfaces_fitted_normalized),[],'all');
        zmin_astigmatism = zmin;
        zmax_astigmatism = zmax;
        wrms_surface_groundtruth = sqrt(sum(zernikeSurfaces_groundtruth_normalized.^2,3));
    else
        zmin = min(zernikeSurfaces_fitted_normalized,[],'all');
        zmax = max(zernikeSurfaces_fitted_normalized,[],'all');
        zmin_astigmatism = zmin;
        zmax_astigmatism = zmax;
    end
elseif strcmp(fitmodel,'xyz-gamma')
    if params.groundtruth_exists
        zmin = min(cat(1,zernikeSurfaces_groundtruth_normalized(:,:,[1 2 4 5 6]),zernikeSurfaces_fitted_normalized(:,:,[1 2 4 5 6])),[],'all');
        zmax = max(cat(1,zernikeSurfaces_groundtruth_normalized(:,:,[1 2 4 5 6]),zernikeSurfaces_fitted_normalized(:,:,[1 2 4 5 6])),[],'all');
        zmin_astigmatism = min(min(cat(1,zernikeSurfaces_groundtruth_normalized(:,:,3),zernikeSurfaces_fitted_normalized(:,:,3)),[],'all'),0);
        zmax_astigmatism = max(cat(1,zernikeSurfaces_groundtruth_normalized(:,:,3),zernikeSurfaces_fitted_normalized(:,:,3)),[],'all');
        wrms_surface_groundtruth = sqrt(sum(zernikeSurfaces_groundtruth_normalized.^2,3));
    else
        zmin = min(zernikeSurfaces_fitted_normalized(:,:,[1 2 4 5 6]),[],'all');
        zmax = max(zernikeSurfaces_fitted_normalized(:,:,[1 2 4 5 6]),[],'all');
        zmin_astigmatism = min(min(zernikeSurfaces_fitted_normalized(:,:,3),[],'all'),0);
        zmax_astigmatism = max(zernikeSurfaces_fitted_normalized(:,:,3),[],'all');
    end
end

wrms_surface_fitted = sqrt(sum(zernikeSurfaces_fitted_normalized.^2,3));

if zmin==zmax
    zmin = zmin - 40;
    zmax = zmax + 40;
end

%% Convert FOV coordinates from [-1,1] to physical coordinates in um.
%Xim = Xn;
%Yim = Yn;

[Xim,Yim] = normalized_to_physical_fov_coordinates(Xn,Yn,params);
[xim,yim] = normalized_to_physical_fov_coordinates(xn,xn,params);

nr_of_zernikes = size(params.aberrations,1);
orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernikes,1);

for jzer = 1:nr_of_zernikes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

%% 
figure('Name','Fitted aberration surfaces')  
for nn = 1:nr_of_zernikes
    hsub = subplot(2,3,nn);
    %contourf(Xim,Yim,zernikeSurfaces_fitted_normalized(:,:,nn));%,zmin:20:zmax)
    imagesc(xim,yim,zernikeSurfaces_fitted_normalized(:,:,nn));
    title(['A(' allxticklabels{nn} ')']);
    %clim([zmin zmax])
    clim([-110 110])
    %clim([-60 60])
    %clim([-50 50])
    %colorbar
    xlabel('x(\mum)')
    ylabel('y(\mum)')
    colormap(redblue)
    colorbar
    ax = gca;
    ax.FontSize = 14;
    axis square
    set(gca,'YDir','normal') 
end


%% Compare to aber map from beads

aber_map_scaled = 1e3*aber_map;
% % downsample
aber_map_downsampled = zeros(25, 25, 23);
for i = 1:23
    originalSlice = aber_map_scaled(:, :, i);
    downsampledSlice = imresize(originalSlice, [25, 25]);
    aber_map_downsampled(:, :, i) = downsampledSlice;
end


aber_map_adjusted = zeros(size(zernikeSurfaces_fitted_normalized));
aber_map_adjusted(:,:,1) = aber_map_downsampled(:,:,1);
aber_map_adjusted(:,:,3) = aber_map_downsampled(:,:,2);
aber_map_adjusted(:,:,4) = aber_map_downsampled(:,:,3);
aber_map_adjusted(:,:,5) = aber_map_downsampled(:,:,4);
aber_map_adjusted(:,:,6) = aber_map_downsampled(:,:,5);
aber_map_adjusted = permute(aber_map_adjusted,[2 1 3]);

%% figure('Name','Aberration surfaces')  
figure('Name','Bead aberration surfaces')  
for nn = 1:nr_of_zernikes
    hsub = subplot(2,3,nn);
    contourf(Xim,Yim,aber_map_adjusted(:,:,nn),zmin:20:zmax)
    clim([zmin zmax])
    ax = gca;
    ax.FontSize = 14;
end


%% Plot 3D surfaces

nr_of_zernikes = size(params.aberrations,1);
orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernikes,1);

for jzer = 1:nr_of_zernikes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

% Zernike surfaces
if params.groundtruth_exists
    zernikeSurfaces_groundtruth = get_zernike_coefficients(Xn,Yn,groundtruth.global,params);
end
zernikeSurfaces_fitted = get_zernike_coefficients(Xn,Yn,theta.global,params);

% Convert aberrations from nm to mlambda
if params.groundtruth_exists
    zernikeSurfaces_groundtruth_normalized = 1e3*zernikeSurfaces_groundtruth/params.lambda;
end
zernikeSurfaces_fitted_normalized = 1e3*zernikeSurfaces_fitted/params.lambda;

zmin = [-200,-200,-200,-200,-200,-200];
zmax = [240,240,240,240,240,240];
view_array = [[-1.8 1 0.8];[1.5 -1 0.4];[1 -1.7 1];
        [-1 -1 0.4];[-1 -1.5 0.4];[2 1 0.8];];

figure('Name','Aberration surfaces')  
for nn = 6%1:nr_of_zernikes
    %hsub = subplot(2,3,nn);    
    %set(hsub, 'Visible', 'off');
    surf(Xim,Yim,aber_map_adjusted(:,:,nn),'FaceAlpha',0.45,'EdgeColor', 'none', 'FaceColor',[0.0745    0.6235    1.0000])
    hold on; axis square;
    surf(Xim,Yim,zernikeSurfaces_fitted_normalized(:,:,nn),'EdgeColor', 'none','FaceAlpha',0.6, 'FaceColor',[1.0000    0.5333         0])
    hold off; axis square;
    %title(['A(' allxticklabels{nn} ')']);
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    zlabel('Zernike coef. (m\lambda)');
    zlim([zmin(nn) zmax(nn)]);
    ax=gca;
    ax.FontSize = 36;
    view(view_array(nn,:));
    ax.LineWidth = 1.5;
    ax.GridAlpha = 0.25;
end

%legend('Fitted surface','Beads','Location','east','Fontsize',20)

%% Compare bead aberrations to single molecule aberrations 
% (Used for figures Lidke data)

clear all
close all

%path_output_data = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\beads\Without_Astigmatism\results_*';
path_output_data = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\beads\With_Astigmatism\results_*';

all_output_files = dir(path_output_data);
nr_of_outfiles = size(all_output_files,1);
fprintf('Found %i files\n',nr_of_outfiles);

theta_combined = [];
roixy_combined = [];
mu_combined = [];
allspots_combined = [];

%for file_i=1:nr_of_outfiles
for file_i=6:9

    filename = all_output_files(file_i).name; 
    foldername = all_output_files(file_i).folder;
    path_output_data_full = fullfile(foldername, filename);    
    
    dataout = load(path_output_data_full,'theta','roixy','mu','allspots','outliers','params');
    
    theta_temp = dataout.theta.local;
    roixy_temp = dataout.roixy;
    mu_temp = dataout.mu;
    allspots_temp = dataout.allspots;
    outliers = dataout.outliers;
    %outliers = [];
    no_outliers = setdiff(1:size(theta_temp,2),outliers);
    params = dataout.params;

    theta_combined = cat(2,theta_combined,theta_temp(:,no_outliers));
    roixy_combined = cat(2,roixy_combined,roixy_temp(:,no_outliers));
    mu_combined = cat(4,mu_combined,mu_temp(:,:,:,no_outliers));
    allspots_combined = cat(4,allspots_combined,allspots_temp(:,:,:,no_outliers));

end

%% Remove spots with large fiterror

errM = get_fiterror(mu_combined,allspots_combined,params);

%idx_err = find(errM(1,:)>50 | errM(2,:)>50 | errM(3,:)>5e4); % Lidke 2D
idx_err = find(errM(1,:)>13 | errM(2,:)>13 | errM(3,:)>5e4); % Lidke 3D
small_err = setdiff(1:size(roixy_combined,2),idx_err);

theta_combined = theta_combined(:,small_err);
roixy_combined = roixy_combined(:,small_err);
mu_combined = mu_combined(:,:,:,small_err);
allspots_combined = allspots_combined(:,:,:,small_err);


%% Plot bead aberrations

pixelsize = params.pixelsize;
xsize = pixelsize*params.imgSizeX;
ysize = pixelsize*params.imgSizeY;

[Xq,Yq] = meshgrid(1:xsize*1e-3);

zernike_coeffs = theta_combined(6:end,:)*1e3/params.lambda;%% Plot aberration surfaces

[~,FOV_coordinates] = get_fov_coordinates(roixy_combined,theta_combined(1,:),theta_combined(2,:),params);
X = FOV_coordinates(1,:);
Y = FOV_coordinates(2,:);

labels = ["A(2,-2)" "A(2,2)" "A(3,-1)" "A(3,1)" "A(4,0)"];
figure
for izer = 4%1:size(params.aberrations,1)
    %hsub = subplot(2,3,izer);
    V = zernike_coeffs(izer,:);
    Vq = griddata(X,Y,V,Xq,Yq,"cubic");

    surf(Xq,Yq,Vq,'FaceAlpha',0.6)
    shading interp
    hold on
    plot3(X,Y,V,"o",'LineWidth',1.5,'MarkerSize',4)
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    zlabel('m\lambda')
    %zlim([-100 100])
    title(labels(izer))
    ax = gca;
    ax.FontSize = 16;
    axis square

end

%% Show beads
n_bead = 2;
dipshow(cat(2,allspots_combined(:,:,:,n_bead),mu_combined(:,:,:,n_bead)),'lin')
colormap parula
diptruesize(1500)
theta_combined(3,n_bead)

%% Plot beads and single molecule aberrations in one plot
%path_aber = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\Data_wo_Astigmatism\fitted_aberrations\archive\rand_init_cppbugfixed_xyz\002_fitted_aberrations_segmentation_thr60_Data0001';
%path_aber = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\Data_wo_Astigmatism\fitted_aberrations\rand_init_cppbugfixed_xyz\002_fitted_aberrations_segmentation_thr60_Data0001.mat';
%path_aber = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\Data_wo_Astigmatism\fitted_aberrations\roi17\fitted_aberrations_0002_segmentation_roi17_offset100_thr30_Data0010.mat';
path_aber = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\Data_w_Astigmatism\fitted_aberrations\001_fitted_aberrations_segmentation_thr30_Data0001.mat';

aber_data = load(path_aber);
%%
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
theta_global = aber_data.theta_full.global;
%theta_global = theta_global.*[-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1]';
params_aber = aber_data.params;
%params_aber = set_parameters_Lidke_3D;
theta_global_init = aber_data.thetainit.global;

x_aber = linspace(-1,1,50);
y_aber = linspace(-1,1,50);
[X_aber,Y_aber] = meshgrid(x_aber,y_aber);

zernikeSurfaces_fitted = get_zernike_coefficients(X_aber,Y_aber,theta_global,params_aber);
zernikeSurfaces_fitted_normalized = 1e3*zernikeSurfaces_fitted/params_aber.lambda;

zernikeSurfaces_init_normalized = 1e3*get_zernike_coefficients(X_aber,Y_aber,theta_global_init,params_aber)/params_aber.lambda;

jzer = [1 0 2 3 4 5];

[Xim,Yim] = normalized_to_physical_fov_coordinates(X_aber,Y_aber,params_aber);

view_array = [[-1.8 1 0.8];[1.5 -1 0.4];[1 -1.7 1];
        [-1 -1 0.4];[-1 -1.5 0.4];[2 1 0.8];];
%%
labels = ["A(2,-2)" "A(2,0)" "A(2,2)" "A(3,-1)" "A(3,1)" "A(4,0)"];
%labels = ["A(2,-2)" "A(2,2)" "A(3,-1)" "A(3,1)" "A(4,0)"];
figure
for izer = [4 5]%[1 2 3 4 5 6]%1:size(params_aber.aberrations,1)
    
    iplot = izer;
    %hsub = subplot(2,3,izer);
    hsub = subplot(1,2,izer-3);
    zs = surf(Xim,Yim,zernikeSurfaces_fitted_normalized(:,:,izer),'FaceAlpha',0.5,'FaceColor','#F80');
    shading interp
    hold on

    if izer~=2
        V = zernike_coeffs(jzer(izer),:);
        Vq = griddata(X,Y,V,Xq,Yq,"cubic");
        surf(Xq,Yq,Vq,'FaceAlpha',0.6)
        shading interp
        hold on
        plot3(X,Y,V,"o",'LineWidth',3,'MarkerSize',6)
        %plot3(X,Y,V,"o",'LineWidth',4,'MarkerSize',8)
        
    end

    set(zs,'Facecolor','#F80');

    xlabel('x (\mum)')
    ylabel('y (\mum)')
    zlabel('Zernike coef. (m\lambda)');
    % 2D + nonastigmatic 3D coeffs
    %zlim([-75 50])

    % 3D astigmatism coeff
    %zlim([-75 150])
    zlim([-50 70])
    title(labels(izer))
    ax=gca;
    ax.FontSize = 24;%32;
    ax.LineWidth = 1.5;
    ax.GridAlpha = 0.25;
    axis square
    view(view_array(izer,:));
end

%% Plot z-location of single molecules in xyz-aberration estimation mode

[Xq,Yq] = meshgrid(1:aber_data.params.imgSizeX*aber_data.params.pixelsize*1e-3);
X = aber_data.roixy(1,:)*aber_data.params.pixelsize*1e-3;
Y = aber_data.roixy(2,:)*aber_data.params.pixelsize*1e-3;
Zq = griddata(X,Y,aber_data.theta_full.local(3,:),Xq,Yq,"cubic");

figure
%int_surf = surf(Xq,Yq,abs(Zq),'FaceAlpha',0.5);
shading interp
hold on
beads = plot3(X,Y,aber_data.theta_full.local(3,:),"o",'LineWidth',2);
xlabel('x (\mum)')
ylabel('y (\mum)')
zlabel('z (nm)')
zlim([-500 800])
title('Sample tilt')
legend('Interpolated surface','bead z-value')
ax = gca;
ax.FontSize = 16;

%% Find fitted aberrations with highest likelihood from multiple random initializations 
% (needs to be with the same spots)

path = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\Data_wo_Astigmatism\fitted_aberrations\roi17\fitted_*';

all_solutions = dir(path);


merit_all = zeros(size(all_solutions,1),1);
merit_all_outliers = zeros(size(all_solutions,1),1);
merit_all_no_outliers = zeros(size(all_solutions,1),1);
nr_of_no_outliers = zeros(size(all_solutions,1),1);

system_aberration_error = zeros(size(all_solutions,1),1);
min_idx = zeros(size(all_solutions,1),1);
theta_global_beads = [0 0 0 0 -2 4 0 0 2 0 -3 -18 13]'; %Lidke 2D

for i=1:size(all_solutions,1)
    i
    path_full = fullfile(all_solutions(i).folder,all_solutions(i).name)
    load(path_full)
    outliers = outliers_full{1}';
    no_outliers = setdiff(1:params.Ncfg,outliers);
    nr_of_no_outliers(i) = numel(no_outliers);
    
    merit_all(i) = sum(meritstore_full.local(:,end))/params.Ncfg;
    merit_all_outliers(i) = sum(meritstore_full.local(outliers,end))/numel(outliers);
    merit_all_no_outliers(i) = sum(meritstore_full.local(no_outliers,end))/numel(no_outliers);

    theta_global_temp = theta_full.global;
    theta_global_diff_1 = theta_global_beads - theta_global_temp;
    theta_global_diff_2 = theta_global_beads - theta_global_temp.*[-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1]';
    sav_1 = get_system_aberration_value(theta_global_diff_1,params);
    sav_2 = get_system_aberration_value(theta_global_diff_2,params);
    [system_aberration_error(i),min_idx(i)] = min([sav_1,sav_2]);

end
%%
figure
x = merit_all;
b = bar(x);

b.FaceColor = 'flat';

ylim([min(x)-0.1*(max(x) - min(x)) max(x)+0.1*(max(x) - min(x))])

xlabel('solution nr')
ylabel('likelihood')

[M,I] = max(x)


%% Find bead indices
xrange = [40,60]*params.imgSizeX/params.pixelsize; %pixels
yrange = [0,20]*params.imgSizeX/params.pixelsize; %pixels
indices = find((roixy_combined(1,:) >= xrange(1) & roixy_combined(1,:) <= xrange(2) & roixy_combined(2,:) >= yrange(1) & roixy_combined(2,:) <= yrange(2) ));

n_bead = indices(4);
dipshow(cat(2,allspots_combined(:,:,:,n_bead),mu_combined(:,:,:,n_bead)),'lin')
%dipshow(allspots_combined(:,:,:,n_bead),'lin')
colormap parula
diptruesize(1500)
theta_combined(3,n_bead)

%% Calculate average astigmatism in central part of FOV

surfaces_center = zernikeSurfaces_fitted_normalized(21:40,21:40,:);
surfaces_center_mean = mean(surfaces_center,[1 2]);
surfaces_center_std = std(surfaces_center,0,[1 2])
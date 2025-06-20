% Caluculate and visualize localtion error with arrows to show the effect
% of field dependent aberrations.
% Difference between fitted values and groundtruth
% Difference between fitting with and without field dependent aberrations

%% Difference between fitted locations and groundtruth
clear all

% % Open locations 1
% path_location_1 = 'read\simulationdata\results\results_fitted_aberrations_1_.mat';
% load(path_location_1,'localizations','outliers','params')
% %Ncfg = size(localizations_full,1);
% %no_outliers = setdiff(1:Ncfg,outliers_full{max_merit_index});
% locations_1 = localizations(:,2:3);

% Open locations 1
%path_location_1 = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\localizations_after_drift_correction\localizations_gain083_3okt.mat';
%path_location_1 = '/home/idroste/Desktop/TUDelft/Data/Lidke/data_analysis/DATA1/Data_wo_Astigmatism/localization/roi9_cppbugfixed/localization_combined.mat';
path_location_1 = '/home/idroste/Desktop/TUDelft/Data/Lidke/data_analysis/DATA1/Data_w_Astigmatism/localization/localization_combined.mat';
data1 = load(path_location_1,'localizations','params');
%locations_1 = localizations(:,2:3);

% % Open locations 2
% path_locations_2 = 'read\simulationdata\input_data\input_spots_1_.mat';
% load(path_locations_2,'groundtruth','roixy')
% [~,locations_2] = get_fov_coordinates(roixy,groundtruth.local(1,:),groundtruth.local(2,:),params);
% locations_2 = 1e3*locations_2';

% % Open locations 2
% path_locations_2 = 'read\simulationdata\results\results_noaber_1_.mat';
% load(path_locations_2,'localizations')
% locations_2 = localizations(:,2:3);

% Open locations 2
%path_locations_2 = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\localizations_after_drift_correction\localizations_gain083_constaber_3okt.mat';
%path_locations_2 = '/home/idroste/Desktop/TUDelft/Data/Lidke/data_analysis/DATA1/Data_wo_Astigmatism/localization/roi9_noaber/localization_combined.mat';
path_locations_2 = '/home/idroste/Desktop/TUDelft/Data/Lidke/data_analysis/DATA1/Data_w_Astigmatism/localization/const_aber/localization_combined.mat';
data2 = load(path_locations_2,'localizations');
%locations_2 = localizations(:,2:3);

%% Find spots that are in both lists
params = data1.params;
localizations_1 = data1.localizations;
localizations_2 = data2.localizations;
[C,index_1,index_2] = intersect(localizations_1(:,1),localizations_2(:,1));

xyz_1 = localizations_1(index_1,2:4);
xyz_2 = localizations_2(index_2,2:4);

%%
% Set parameters
imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;
pixelsize = params.pixelsize;

location_difference = xyz_1 - xyz_2;

%%

% define FOV patches
Npatch_xy = 15;
Npatch_z = 1;

zmin = min(cat(1,localizations_1(:,4),localizations_2(:,4)),[],1);
zmax = max(cat(1,localizations_1(:,4),localizations_2(:,4)),[],1);

minxyz = [0 0 -1000];
maxxyz = [imgSizeX*pixelsize imgSizeY*pixelsize 1000];
Xpatch = minxyz(1)+(0:Npatch_xy)*(maxxyz(1)-minxyz(1))/Npatch_xy;
Ypatch = minxyz(2)+(0:Npatch_xy)*(maxxyz(2)-minxyz(2))/Npatch_xy;
Zpatch = minxyz(3)+(0:Npatch_z)*(maxxyz(3)-minxyz(3))/Npatch_z;

mean_location_difference = zeros(3,Npatch_xy,Npatch_xy,Npatch_z);
spots_per_patch = zeros(Npatch_xy,Npatch_xy,Npatch_z);

for ii = 1:Npatch_xy
    ii
    for jj = 1:Npatch_xy
        for kk=1:Npatch_z
            FOVfilter = (xyz_1(:,1)>Xpatch(ii))&(xyz_1(:,1)<Xpatch(ii+1))&...
                        (xyz_1(:,2)>Ypatch(jj))&(xyz_1(:,2)<Ypatch(jj+1))&...
                        (xyz_1(:,3)>Zpatch(kk))&(xyz_1(:,3)<Zpatch(kk+1));
    
            mean_location_difference(:,ii,jj,kk) = mean(location_difference(FOVfilter,:));
            spots_per_patch(ii,jj,kk) = size(location_difference(FOVfilter,:),1);
        end
    end
end

%% Plot arrows (xy, 1 img for all z)

scalefac = 6; %0.2;
x_legend = -5;
y_legend = -4;
u_legend = 2*scalefac;
v_legend = 0;

figure
x = 1e-3*Xpatch(1:end-1);
y = 1e-3*Ypatch(1:end-1);
[xx,yy] = meshgrid(x,y);
mean_location_difference_scaled = -1*scalefac*mean_location_difference;
%quiver3(xx(:)',yy(:)',zeros(size(xx(:)')),mean_location_difference_scaled(1,:),mean_location_difference_scaled(2,:),mean_location_difference_scaled(3,:),0,'LineWidth',2)
quiver(x_legend,y_legend,u_legend,v_legend,0,'LineWidth',2.0,'color',"#0072BD","MaxHeadSize",10)
hold on
quiver(xx(:)',yy(:)',mean_location_difference_scaled(2,:),mean_location_difference_scaled(1,:),0,'LineWidth',2,'AutoScale','off','color',"#0072BD","MaxHeadSize",10)
axis square
xlim([-6 98])
ylim([-6 98])
xlabel('x (\mum)')
ylabel('y (\mum)')
ax=gca;
ax.FontSize = 16;

% figure
% x = 1e-3*Xpatch(1:end-1);
% y = 1e-3*Ypatch(1:end-1);
% [xx,yy] = meshgrid(x,y);
% mean_location_difference_z_scaled = 0.4*mean_location_difference;
% quiver(xx(:)',yy(:)',mean_location_difference_z_scaled(1,:),mean_location_difference_z_scaled(2,:),0,'LineWidth',2)
% xlabel('x (\mum)')
% ylabel('y (\mum)')

%% xyz-location difference, for all z in one image

scalefac  = 2e-1;

%x = 1e-3*Xpatch(1:end-1);
%y = 1e-3*Ypatch(1:end-1);
x = linspace(1,1e-3*params.pixelsize*params.imgSizeX,Npatch_xy);
y = linspace(1,1e-3*params.pixelsize*params.imgSizeY,Npatch_xy);
z = zeros(Npatch_xy,Npatch_xy);
%x_legend = -30;
%y_legend = -15;
%u_legend = 50*scalefac;
%v_legend = 0;

figure
U = scalefac*squeeze(mean_location_difference(2,:,:,:));
V = scalefac*squeeze(mean_location_difference(1,:,:,:));
W = squeeze(mean_location_difference(3,:,:,:));
%quiver(x_legend,y_legend,u_legend,v_legend,0,'LineWidth',1.5,'color',"#0072BD","MaxHeadSize",15)
%hold on
quiver3(x,y,z,U,V,W,0,'LineWidth',1.5,'color',"#0072BD","MaxHeadSize",0.05)
%xlim([-100 1100])
%ylim([-100 1100])
zlim([-10 1e-3*params.pixelsize*params.imgSizeX-10]) % For correct scaling between xy and z.
xlabel('x (\mum)')
ylabel('y (\mum)')
zlabel('z (nm)')
%text(-30,-25,'50 nm','FontSize',8)
%text(50,-20,strcat('Number of spots: '," ",num2str(sum(spots_per_patch(:,:,iz),"all"))),'FontSize',10)
ax=gca;
ax.FontSize = 16;



%% Make movie of xy at different z-heights

moviedir = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\location_arrows\';
writerObjPSF = VideoWriter(strcat(moviedir,'location_fitted_vs_const_aberrations','.avi'));
writerObjPSF.FrameRate = 1;
open(writerObjPSF);

Zpatchmid = (Zpatch(1:end-1)+Zpatch(2:end))/2;
scalefac = 0.25;

x = 1e-3*Xpatch(1:end-1);
y = 1e-3*Ypatch(1:end-1);
x_legend = -30;
y_legend = -15;
u_legend = 50*scalefac;
v_legend = 0;

for iz = 1:numel(Zpatchmid)
    figure
    U = squeeze(mean_location_difference_z(1,:,:,iz));
    V = squeeze(mean_location_difference_z(2,:,:,iz));
    U = scalefac*U;
    V = scalefac*V;
    quiver(x_legend,y_legend,u_legend,v_legend,0,'LineWidth',1.5,'color',"#0072BD","MaxHeadSize",15)
    hold on
    quiver(x,y,U,V,0,'LineWidth',1.5,'color',"#0072BD","MaxHeadSize",0.2)
    xlim([-40 190])
    ylim([-40 190])
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    title(strcat('z = ',num2str(Zpatchmid(iz)),' nm'),'FontSize',16)
    text(-30,-25,'50 nm','FontSize',8)
    %text(50,-20,strcat('Number of spots: '," ",num2str(sum(spots_per_patch(:,:,iz),"all"))),'FontSize',10)
    ax=gca;
    ax.FontSize = 16;

    frame = getframe(gcf);
    writeVideo(writerObjPSF,frame);
end
close(writerObjPSF);

%%  Make movie of xyz at different z-heights

moviedir = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\location_arrows\';
writerObjPSF = VideoWriter(strcat(moviedir,'location_fitted_vs_const_aberrations_3D','.avi'));
writerObjPSF.FrameRate = 1;
open(writerObjPSF);

Zpatchmid = (Zpatch(1:end-1)+Zpatch(2:end))/2;
scalefac = 1;

%x = 1e-3*Xpatch(1:end-1);
%y = 1e-3*Ypatch(1:end-1);
x = linspace(1,1000,Npatch_xy);
y = linspace(1,1000,Npatch_xy);
z = zeros(Npatch_xy,Npatch_xy);
%x_legend = -30;
%y_legend = -15;
%u_legend = 50*scalefac;
%v_legend = 0;


for iz = 1:numel(Zpatchmid)
    figure
    U = squeeze(mean_location_difference_z(1,:,:,iz));
    V = squeeze(mean_location_difference_z(2,:,:,iz));
    W = squeeze(mean_location_difference_z(3,:,:,iz));
    %quiver(x_legend,y_legend,u_legend,v_legend,0,'LineWidth',1.5,'color',"#0072BD","MaxHeadSize",15)
    %hold on
    quiver3(x,y,z,U,V,W,0,'LineWidth',1.5,'color',"#0072BD","MaxHeadSize",0.2)
    xlim([-100 1100])
    ylim([-100 1100])
    zlim([-600 600])
    xlabel('x (nm)')
    ylabel('y (nm)')
    zlabel('z (nm)')
    title(strcat('z = ',num2str(Zpatchmid(iz)),' nm'),'FontSize',16)
    %text(-30,-25,'50 nm','FontSize',8)
    %text(50,-20,strcat('Number of spots: '," ",num2str(sum(spots_per_patch(:,:,iz),"all"))),'FontSize',10)
    ax=gca;
    ax.FontSize = 16;

    frame = getframe(gcf);
    writeVideo(writerObjPSF,frame);
end
close(writerObjPSF);

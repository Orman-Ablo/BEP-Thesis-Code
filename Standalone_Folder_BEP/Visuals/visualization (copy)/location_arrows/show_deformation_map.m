% Show local deformation when not fitting aberrations compared to fitting
% aberrations.

%% Difference between fitted locations and groundtruth
clear all

% Open locations 1
path_location_1 = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\localizations_after_drift_correction\localizations_gain083_3okt.mat';
data1 = load(path_location_1,'localizations','params');
%locations_1 = localizations(:,2:3);

% Open locations 2
path_locations_2 = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\localizations_after_drift_correction\localizations_gain083_constaber_3okt.mat';
data2 = load(path_locations_2,'localizations');
%locations_2 = localizations(:,2:3);

%% Find spots that are in both lists
params = data1.params;
localizations_1 = data1.localizations;
localizations_2 = data2.localizations;
[C,index_1,index_2] = intersect(localizations_1(:,1),localizations_2(:,1));

xyz_1 = localizations_1(index_1,2:4);
xyz_2 = localizations_2(index_2,2:4);

location_difference = xyz_1 - xyz_2;

%% Compute average per patch
imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;
pixelsize = params.pixelsize;

Npatch_xy = 25;
Npatch_z = 25;

zmin = min(cat(1,localizations_1(:,4),localizations_2(:,4)),[],1);
zmax = max(cat(1,localizations_1(:,4),localizations_2(:,4)),[],1);

minxyz = [0 0 -400];
maxxyz = [imgSizeX*pixelsize imgSizeY*pixelsize 400];
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

%% Calculate deformation
% Deformation maps dx/dx, dx/dy, dx/dz etc. (9 in total)
% Derivative: f'(x) = (f(x+h) - f(x))/h 
% Set parameters

derivative_x = (mean_location_difference(:,2:Npatch_xy,:,:) - mean_location_difference(:,1:Npatch_xy-1,:,:))/(imgSizeX*pixelsize/Npatch_xy);
derivative_y = (mean_location_difference(:,:,2:Npatch_xy,:) - mean_location_difference(:,:,1:Npatch_xy-1,:))/(imgSizeX*pixelsize/Npatch_xy);
derivative_z = (mean_location_difference(:,:,:,2:Npatch_z) - mean_location_difference(:,:,:,1:Npatch_z-1))/((maxxyz(3) - minxyz(3))/Npatch_z);

% % Replace NaN (no locs in this patch) by zero
% nanIdx_x = isnan(derivative_x);
% nanIdx_y = isnan(derivative_y);
% nanIdx_z = isnan(derivative_z);
% 
% derivative_x(nanIdx_x) = 0;
% derivative_y(nanIdx_y) = 0;
% derivative_z(nanIdx_z) = 0;

%% Show deformation maps
imagesc(squeeze(derivative_z(3,:,:,5)))
cb = colorbar;
ylabel(cb,'nm','Rotation',270,'FontSize',16)
ax=gca;
ax.FontSize = 16;

%% Make a movie of deformation maps

moviedir = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\location_arrows\';
writerObjPSF = VideoWriter(strcat(moviedir,'deformation_map_dxdz','.avi'));
writerObjPSF.FrameRate = 1;
open(writerObjPSF);

Zpatchmid = (Zpatch(1:end-1)+Zpatch(2:end))/2;
scalefac = 1;

%x = 1e-3*Xpatch(1:end-1);
%y = 1e-3*Ypatch(1:end-1);
%x = linspace(1,1000,Npatch_xy);
%y = linspace(1,1000,Npatch_xy);
%z = zeros(Npatch_xy,Npatch_xy);
%x_legend = -30;
%y_legend = -15;
%u_legend = 50*scalefac;
%v_legend = 0;

xrange = [0 imgSizeX*pixelsize*1e-3];
yrange = [0 imgSizeY*pixelsize*1e-3];


for iz = 1:Npatch_z-1
    imagesc(xrange,yrange,squeeze(derivative_z(1,:,:,iz)))
    clim([-3 3]);
    colormap([1,1,1;hsv]);
    cb = colorbar;
    cb.Ticks = -3:3;
    ylabel(cb,'nm','Rotation',270,'FontSize',16)
    ax=gca;
    ax.FontSize = 16;
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    title(["Deformation dx/dz",strcat('z = ',num2str(Zpatchmid(iz)),' nm')]);
    %text(-30,-25,'50 nm','FontSize',8)
    %text(50,-20,strcat('Number of spots: '," ",num2str(sum(spots_per_patch(:,:,iz),"all"))),'FontSize',10)
    ax=gca;
    ax.FontSize = 14;

    frame = getframe(gcf);
    writeVideo(writerObjPSF,frame);
end
close(writerObjPSF);
%%

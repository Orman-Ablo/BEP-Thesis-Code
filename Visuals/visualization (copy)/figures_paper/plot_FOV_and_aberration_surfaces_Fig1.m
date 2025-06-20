% Plot aberration surfaces above FOV with spots
close all
clear all

params = set_parameters_xy_NAT_Fig1;

Ncfg_total = 100;
N_level = 2500;
b_level = 10;
xy_range = 2;
z_range = [-500,500];% in nm
pixelsize = params.pixelsize;

params.add_newnoise = false;

% Change aberrations in generate_NAT_data if you do not want random
% aberrations
[groundtruth,roixy,allspots,dmudthetastore,fov_coordinates,zernikeCoefficients,Wrms,framelist,ID] = generate_NAT_data(Ncfg_total,xy_range,z_range,N_level,b_level,params,0,0);

%dipshow(log(allspots),'lin'); colormap parula; diptruesize(1000);
dipshow(poissrnd(allspots),'lin'); colormap parula; diptruesize(3000);

xn = linspace(-1,1,50);
yn = linspace(-1,1,50); 
[Xn,Yn] = meshgrid(xn,yn);
zernikeSurfaces = get_zernike_coefficients(Xn,Yn,groundtruth.global, params);

imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;
imgSizeZ = params.Mz;
%full_image = zeros(imgSizeY,imgSizeX,imgSizeZ);
full_image = 45*ones(imgSizeY,imgSizeX,imgSizeZ);
Mx = params.Mx;
My = params.My;
Mz = params.Mz;

xmax = pixelsize*imgSizeX*1e-03;
xmin = 0;
ymax = pixelsize*imgSizeY*1e-03;
ymin = 0;

for jcfg=1:Ncfg_total
    full_image(roixy(1,jcfg):roixy(1,jcfg)+Mx-1, roixy(2,jcfg):roixy(2,jcfg)+My-1,:) = full_image(roixy(1,jcfg):roixy(1,jcfg)+Mx-1, roixy(2,jcfg):roixy(2,jcfg)+My-1,:) + allspots(:,:,:,jcfg);  
end

full_image = poissrnd(2*full_image);

figure
imagesc(log(max(full_image,0)),[4.5 5.4])
axis square

%%
[Xim,Yim] = normalized_to_physical_fov_coordinates(Xn,Yn,params);
nr_of_zernikes = size(params.aberrations,1);
orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernikes,1);

for jzer = 1:nr_of_zernikes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

[fov_coordinates_normalized,fov_coordinates_physical] = get_fov_coordinates(roixy,groundtruth.local(1,:),groundtruth.local(2,:),params);

for nn = 1:nr_of_zernikes
        figure1 = figure;
        axes1 = axes('Parent',figure1);
        %hold(axes1,'on');
        %subplot(2,3,nn)    
        surf(1*Xim,1*Yim,zernikeSurfaces(:,:,nn),'FaceAlpha',0.75);
        hold on; axis square; shading interp;
        % Plot FOV in plane z=0:
        %surf([0 1*xmax], [0 1*ymax],repmat(min(zernikeSurfaces(:,:,nn),[],'all'),[2 2]), 15*log(full_image)-55,'facecolor','texture');
        scatter3(fov_coordinates_physical(1,:),fov_coordinates_physical(2,:),zernikeCoefficients(:,nn),100,'o','filled','MarkerFaceColor',[0 0 0]);
        %title(['Anm(' allxticklabels{nn} ')'],FontSize=20);
        zlim([min(zernikeSurfaces(:,:,nn),[],'all') max(zernikeSurfaces(:,:,nn),[],'all')+1])
        xlabel('x (\mum)','FontSize',20);
        ylabel('y (\mum)','FontSize',20);
        title(['A(' allxticklabels{nn} ')']);
        %zl = zlabel(["A_2^2", "(m\lambda)"],'Rotation',0);
        %zl.Position(2) = zl.Position(2) + 10; 
        ax=gca;
        ax.FontSize = 32;
        ax.GridAlpha = 0.35; % maximum line opacity is 1
        view([-67 34.5]);
 end
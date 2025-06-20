% This script generates synthetic spots with aberrations according to Nodal
% Aberration Theory (NAT). Then noise is added and the spots are
% saved and plotted. Four spots are located in the 4 outer corners of the
% FOV, the other ones are randomly distributed over the FOV.
% It is also possible to only plot previously generated data.
%
% Autor: Isabel Droste, TU Delft, 2022

close all;
%clearvars -except allspots
%addpath(genpath('C:\Users\idroste\Desktop\TUDelft\Code\vecfitcpu_vortex'));

%% Settings
params = set_parameters_xy_NAT;
generate_data = true;
save_data = false;
plot_data = true;

imgSizeX = params.imgSizeX; % in nr. of pixels
imgSizeY = params.imgSizeY;
pixelsize = params.pixelsize;
Mx = params.Mx;
My = params.My;
Mz = params.Mz;

% Local parameter settings
xy_range = 2; % Size of area where location is randomly generated, in nr. of pixels.
N_level = 2500;
b_level = 0;

% Global parameter settings
%gamma.defocus = 
%gamma.primary_astigmatism = [-0.49 4.46 0.022 0.015 2.2e-04];
%gamma.primary_coma = [18.4 2.8 -0.06 0];
%gamma.primary_sperical = -19.6;

% For poster PSFs, zero aberrations at center
xmax = pixelsize*imgSizeX*1e-03;
xmin = 0;
xcenter = (xmax - xmin)/2;
ymax = pixelsize*imgSizeY*1e-03;
ymin = 0;
ycenter = (ymax - ymin)/2;

astigmatism_b = 6;
astigmatism_c = 0;
astigmatism_d = 0;
%astigmatism_b = 0.0;
%astigmatism_c = 0.0;
%astigmatism_d = 0.1;
astigmatism_a = -astigmatism_b*xcenter - astigmatism_c*ycenter - astigmatism_d*2*xcenter*ycenter;
astigmatism_a_prime = astigmatism_c*xcenter - astigmatism_b*ycenter - astigmatism_d*(ycenter^2-xcenter^2);

gamma.primary_astigmatism = [astigmatism_a astigmatism_a_prime astigmatism_b astigmatism_c astigmatism_d];
%gamma.primary_astigmatism = [0 0 0 0 0];

% calculate e and e' such that coma 0 in the center of the FOV
coma_coefficient_f = 0.6;
coma_coefficient_g = 0.6;
coma_coefficient_f = 2.4;
coma_coefficient_g = 2.4;
coma_coefficient_e = -coma_coefficient_f * xcenter - coma_coefficient_g * ycenter;
coma_coefficient_e_prime = coma_coefficient_f * ycenter - coma_coefficient_g * xcenter;
gamma.primary_coma = [coma_coefficient_e coma_coefficient_e_prime coma_coefficient_f coma_coefficient_g];
%gamma.primary_coma = [0 0 0 0];

gamma.primary_sperical = 0;

source_file = ''; % Only needed when generate_data = false.
dest_file = 'N:\tnw\IST\QI\users\idroste\Data\synthetic_data\test_beads\bead_data.mat'; % Where data will be saved.

% Include dipimage
addpath('C:\Users\idroste\Desktop\TUDelft\Code\diplib\share\DIPimage');
setenv('PATH',['C:\Users\idroste\Desktop\TUDelft\Code\diplib\bin',';',getenv('PATH')]);

%% Generate data

    
    %[groundtruth,roixy,allspots,global_coordinates,zernikeCoefficients,Wrms] = generate_NAT_data(Ncfg,xy_range,N_level,b_level,gamma,params);
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%imgsize

%% Generate local and global parameters
Ncfg = 100;
if generate_data
    fprintf('Generating corner spots\n')
    % Generate ROIs in corner points 
    roixy = zeros(Ncfg,2);

    % Create 4 corner spots
    xmin_pix = 1;
    ymin_pix = 1;
    xmax_pix = imgSizeX - Mx;
    ymax_pix = imgSizeY - My;

    roixy(1,1) = xmin_pix;
    roixy(1,2) = ymin_pix;

    roixy(2,1) = xmin_pix;
    roixy(2,2) = ymax_pix;

    roixy(3,1) = xmax_pix;
    roixy(3,2) = ymin_pix;

    roixy(4,1) = xmax_pix;
    roixy(4,2) = ymax_pix;

    % Create randomly distributed spots
    roixy(5:end,1) = randi(xmax_pix,Ncfg-4,1); % pixel units
    roixy(5:end,2) = randi(ymax_pix,Ncfg-4,1);
    
    % Generate x,y uniformly in the range [-pixelsize, pixelsize]
    groundtruth_x = xy_range*pixelsize*(rand(Ncfg,1)-1/2);
    groundtruth_y = xy_range*pixelsize*(rand(Ncfg,1)-1/2);
    groundtruth_N = N_level*ones(Ncfg,1);
    groundtruth_b = b_level*ones(Ncfg,1);
    
    groundtruth.local = [groundtruth_x,groundtruth_y,groundtruth_N,groundtruth_b];
    gamma = zeros(14,1);
    groundtruth.global = gamma;
    
    %% Generate spots
    
    % Transform ROI coordinates (nm) to global coordinates (um)
    global_coordinates_x = (roixy(:,1)*pixelsize + groundtruth_x)*1e-03;
    global_coordinates_y = (roixy(:,2)*pixelsize + groundtruth_y)*1e-03;
    global_coordinates = cat(2,global_coordinates_x,global_coordinates_y);
    
    % Combine into one array
    %global_coordinates_xy = zeros(..);
    
    % Calculate aberrations
    zernikeCoefficients = get_zernike_coefficients(global_coordinates_x,global_coordinates_y,gamma, params);
    
    % Generate PSFs
    allspots = zeros(Mx,My,Mz,Ncfg);
    for jcfg=1:Ncfg
        params.xemit = groundtruth_x(jcfg);
        params.yemit = groundtruth_y(jcfg);
        [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = ...
        get_pupil_matrix(params,zernikeCoefficients(jcfg,:));
        [FieldMatrix,FieldMatrixDerivatives] = ...
        get_field_matrix_derivatives(params,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
        [PSF,~] = get_psfs_derivatives(params,PupilMatrix,FieldMatrix,FieldMatrixDerivatives);
        PSF = groundtruth_N(jcfg)*PSF+groundtruth_b(jcfg);
        allspots(:,:,:,jcfg) = PSF;
    end
    
    % Compute Wrms value. It is equal to the sum of the squared zernike
    % coefficients.
    Wrms = sqrt(sum(zernikeCoefficients.^2,2));
    Wrms = Wrms/params.lambda; % convert to lambda units.

    % Add poisson noise.
    %allspots = poissrnd(allspots);
    
    % ! use imgSize instead
    %params.FOV = imgsizeX; 

else
    loaded_data = load(source_file);
    groundtruth = loaded_data.groundtruth;
    roixy = loaded_data.roixy;
    allspots = loaded_data.allspots;
end

% Create image with all spots in 1 image.
imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;
imgSizeZ = params.Mz;
full_image = zeros(imgSizeY,imgSizeX,imgSizeZ);
Mx = params.Mx;
My = params.My;
Mz = params.Mz;

for jcfg=1:Ncfg
    % roixy is lowest coordinate pair!
    %full_image(roixy(jcfg,1):roixy(jcfg,1)+Mx-1, roixy(jcfg,2):roixy(jcfg,2)+My-1) = allspots(:,:,:,jcfg); 
    full_image(roixy(jcfg,1):roixy(jcfg,1)+Mx-1, roixy(jcfg,2):roixy(jcfg,2)+My-1,:) = full_image(roixy(jcfg,1):roixy(jcfg,1)+Mx-1, roixy(jcfg,2):roixy(jcfg,2)+My-1,:) + allspots(:,:,:,jcfg); 
end

% Add noise
full_image = poissrnd(full_image);
%disp(size(full_image));

% Calculate global image coordinates -> get it from generate_NAT_data


%% Save/plot data.
if save_data
    save(dest_file,'allspots','roixy','framelist','params','groundtruth','zernikeCoefficients','Wrms');
end

if plot_data

    % Plot total FOV
    %dipshow(full_image,'log');
    %diptruesize(60)
    xmax = params.imgSizeX*params.pixelsize*1e-3;
    ymax = params.imgSizeY*params.pixelsize*1e-3;
    surf([0 xmax], [0 ymax],repmat(0,[2 2]), log(full_image),'facecolor','texture');
    xl = xlabel('x (\mum)','FontSize',32);
    %xl.Position(3) = xl.Position(3) - 0.1; 
    yl = ylabel('y (\mum)','FontSize',32);
    %yl.Position(2) = yl.Position(2) - 3; 
    set(gca,'ZTick',[]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    zlim([0 1])
    ax=gca;
    ax.FontSize = 32;
    view([-67 34.5]);
    colormap parula

    if true
    % Plot the 4 corner spots in separate images
    for i=1:4
        dipshow(allspots(:,:,:,i),'log');
        diptruesize(2000)
        colormap parula
    end

    
    % Plot all Zernike surfaces
    xim = (params.imgSizeX*params.pixelsize*1e-3)*linspace(0,1,50);
    yim = (params.imgSizeY*params.pixelsize*1e-3)*linspace(0,1,50);
    [Xim,Yim] = meshgrid(xim,yim);

    zernikeSurfaces = get_zernike_surfaces(Xim,Yim,gamma,params);
    nr_of_zernikes = size(params.aberrations,1);

    orders = params.aberrations(:,1:2);
    allxticklabels = cell(nr_of_zernikes,1);
    for jzer = 1:nr_of_zernikes
        allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
    end

    figure
        set(gcf,'Position',[25 90 1255 880])
    for nn = 1:nr_of_zernikes
        subplot(2,3,nn)    
        surf(Xim,Yim,zernikeSurfaces(:,:,nn),'FaceAlpha',0.75)
        hold on; axis square; shading interp;
        title(['Anm(' allxticklabels{nn} ')']);
        xlabel('x (\mum)');
        ylabel('y (\mum)');
        zlabel('Aberration coefficient (m\lambda)');
    %     xlim([-60 60]);
    %     ylim([-60 60]);
    %     xticks(-60:30:60)
    %     yticks(-60:30:60)
    %    xlim([ -110 110]);
    %    ylim([-110 110]);
    %    xticks(-110:55:110)
    %    yticks(-110:55:110)            
    end

    % Plot coma Zernike surfaces in separate figures
        
    for nn = 1:2
        figure1 = figure;
        axes1 = axes('Parent',figure1);
        %hold(axes1,'on');
        %set(gcf,'Position',[25 90 1255 880])
        %subplot(2,3,nn)    
        surf(4*Xim,4*Yim,zernikeSurfaces(:,:,nn),'FaceAlpha',0.8)
        hold on; axis square; shading interp;
        surf([0 4*xmax], [0 4*ymax],repmat(min(-60,[],'all'),[2 2]), 8.5*log(full_image),'facecolor','texture');
        %scatter3(global_coordinates(:,1),global_coordinates(:,2),zernikeCoefficients(:,nn),25,'o','filled','MarkerFaceColor',[0 0 0]);
        %title(['Anm(' allxticklabels{nn} ')'],FontSize=20);
        xlabel('x (\mum)','FontSize',20);
        ylabel('y (\mum)','FontSize',20);
        %zl = zlabel(['Anm(' allxticklabels{nn} ') (m\lambda)']);
        zl = zlabel(["A_2^-^2", "(m\lambda)"],'Rotation',0);
        zl.Position(2) = zl.Position(2) + 10; 
        ax=gca;
        ax.FontSize = 32;
        ax.GridAlpha = 0.35; % maximum line opacity is 1
        view([-67 34.5]);
    end

    end
    
end




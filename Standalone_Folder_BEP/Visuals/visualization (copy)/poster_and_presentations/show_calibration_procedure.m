% This script plots a fitted zernike surface through fitted bead
% aberrations.
close all;
%clearvars;

% load bead fits
base_directory = 'N:\tnw\IST\QI\users\idroste\Data\data_bead_fitting_22_02_15\processed\zernike12_50it\';
%base_directory = 'N:\tnw\IST\QI\users\idroste\Data\fitted_data_qingru\';
addpath(base_directory);
%base_directory = 'C:\Users\idroste\Desktop\TUDelft\Data\data_bead_fitting_22_02_15\fitted_data_qingru\';

%new
%z1=load(strcat(base_directory,'new_region1_141f_0214.mat'));
%z2=load(strcat(base_directory,'new_region2_141f_0214.mat'));
%z3=load(strcat(base_directory,'new_region3_141f_0214.mat'));
%z4=load(strcat(base_directory,'new_region4_141f_0214.mat'));
%z5=load(strcat(base_directory,'new_region5_141f_0214.mat'));
%z6=load(strcat(base_directory,'new_region6_141f_0215.mat'));
%z7=load(strcat(base_directory,'new_region7_141f_0215.mat'));
%z8=load(strcat(base_directory,'new_region8_141f_0215.mat'));
%z9=load(strcat(base_directory,'new_region9_141f_0215.mat'));
%z10=load(strcat(base_directory,'new_region10_141f_0215.mat'));


%% old
%z1 = load(strcat(base_directory,'localization_old_region_1.mat'));
%z = load(strcat(base_directory,'localization_old_region_5.mat'));

%z1=load(strcat(base_directory,'old_region1_141f_0214.mat'));
%z2=load(strcat(base_directory,'old_region2_141f_0214.mat'));
%z3=load(strcat(base_directory,'old_region3_141f_0214.mat'));
%z4=load(strcat(base_directory,'old_region4_141f_0214.mat'));
%z5=load(strcat(base_directory,'old_region5_141f_0214.mat'));
%z6=load(strcat(base_directory,'old_region6_141f_0215.mat'));
%z7=load(strcat(base_directory,'old_region7_141f_0215.mat'));
%z8=load(strcat(base_directory,'old_region8_141f_0215.mat'));
%z9=load(strcat(base_directory,'old_region9_141f_0215.mat'));
%z10=load(strcat(base_directory,'old_region10_141f_0215.mat'));

%theta=[z1.theta z2.theta z3.theta z4.theta z5.theta z6.theta z7.theta z8.theta z9.theta z10.theta];
%theta=[z1.theta z5.theta];
theta = z.theta;

params=z.params;

addpath('funnat');
%%
%params = set_parameters_zstack_bead;
% params.flg_parallel = 1;
% read data (non-segmented beda)
% file = '\bead\demo_region1_61f';
% datain = ['read' file '.tif'];
% [allspots,roixy,framelist,params] = get_segmentation(params,datain,200);
% xn = params.pixelsize/1E3*(roixy(:,1)-params.FOV/2);
% yn = params.pixelsize/1E3*(roixy(:,2)-params.FOV/2);
%%

a=size(theta);
params.Ncfg=a(:,2);
fov = params.FOV;
pixelsize = params.pixelsize;

% estimated zernike amplitudes
zers = theta(6:end,:)';

% generate grid

% xn = pixelsize/1E3*theta(1,:)'; % (um)
% yn = pixelsize/1E3*theta(2,:)'; % (um)
% roixy(:,1)= roixy(:,1)-params.FOV/2;
% roixy(:,2)= roixy(:,2)-params.FOV/2;
% theta = roi2fov(theta,roixy,params);

xn = pixelsize/1E3*theta(1,:)';
yn = pixelsize/1E3*theta(2,:)';

figure
set(gca,'FontSize',32);
    hold on
    plot(xn,yn,'k+')
      xlim([-70  70]);
     ylim([-70 70]);
     xticks(-70:35:70)
    yticks(-70:35:70)
   xlabel('x (\mum)');
    ylabel('y (\mum)');
    
xim = pixelsize/1E3*((1:fov)-fov/2);
yim = pixelsize/1E3*((1:fov)-fov/2);
xim = xim(1:10:end); % downsample grid
yim = yim(1:10:end);
[Xim,Yim] = meshgrid(xim,yim);

% get nat coefficients, surfaces and predictions
[RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6] = get_natCoefficients_linear(xn,yn,zers,params);
[Zsurface] = get_natSurfaces_linear(Xim,Yim,RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6);
[Zpredict] = get_natPredictions_linear(xn,yn,RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6);

%% Plot
mean(theta(6:17,:),2)
numzers = params.numparams-5-3;
orders = params.aberrations(:,1:2);
allxticks = 1:numzers;
allxticklabels = cell(numzers,1);
for jzer = 1:numzers
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

% plot fitted zernike surfaces and data points (black dots)
figure
    set(gcf,'Position',[25 90 1255 880])
for nn = 1:numzers
    %size(zers,2)
    subplot(3,4,nn)
    surf(xim,yim,Zsurface(:,:,nn),'FaceAlpha',0.75)
     hold on; axis square; shading interp;
    scatter3(xn,yn,zers(:,nn),2,'o','filled','MarkerFaceColor',[0 0 0])
    [r2,rmse] = rsquare(zers(:,nn),Zpredict(:,nn));
    title(strcat(['Znm(' allxticklabels{nn} '), RMSE = ' num2str(rmse,2)]));
    xlabel('x (\mu m)');
    ylabel('y (\mu m)');
    zlabel('Aberration coefficient (m\lambda)');
     xlim([-60 60]);
     ylim([-60 60]);
     xticks(-60:30:60)
     yticks(-60:30:60)
    xlim([ -110 110]);
    ylim([-110 110]);
    xticks(-110:55:110)
    yticks(-110:55:110)
         
end

% Plot only the beads in 3D for a single Zernike mode
%{
nn = 2;
figure
set(gcf,'Position',[25 90 1255 880])
hold on; axis square; shading interp;
scatter3(xn,yn,zers(:,nn), 6,'o','filled','MarkerFaceColor',[0 0 0])
title(['Znm(' allxticklabels{nn} ')']);
xlabel('x (\mu m)');
ylabel('y (\mu m)');
zlabel('Aberration coefficient (m\lambda)');
xlim([-60 60]);
ylim([-60 60]);
xticks(-60:30:60)
yticks(-60:30:60)
xlim([ -110 110]);
ylim([-110 110]);
xticks(-110:55:110)
yticks(-110:55:110)
%}

% Plot the beads + surface for a single Zernike mode
% For each bead add an error bar of size equal to the CRLB value of that
% bead.

for nn = 3:4
    figure
    %set(gcf,'Position',[25 90 1255 880])
    surf(xim,yim,Zsurface(:,:,nn),'FaceAlpha',0.75)
    hold on; axis square; shading interp;
    
    % Plot a dot for each bead.
    scatter3(xn,yn,zers(:,nn), 40,'black','o','filled')
    %scatter3(xn,yn,Zpredict(:,nn), 80,'black','o','filled')

    xlabel('x (\mu m)');
    ylabel('y (\mu m)');
    zlabel('Aberration coefficient (m\lambda)');
    zl = zlabel(["A_3^-^1", "(m\lambda)"],'Rotation',0);
    zl.Position(2) = zl.Position(2) + 40; 
    xlim([ -110 110]);
    ylim([-110 110]);
    zlim([-5 25]);
    xticks(-110:55:110)
    yticks(-110:55:110)
    ax=gca;
    ax.GridLineStyle = '-';
    ax.GridColor = 'k';
    ax.GridAlpha = 0.35; % maximum line opacity is 1
    ax.FontSize = 32;
end

% Check if the error between the measured and predicted aberration values
% is normally distributed by making a Q-Q plot.

%{
for nn = 1:12
    figure
    error_measured_predict = zers(:,nn) - Zpredict(:,nn);
    error_measured_predict_normalized = error_measured_predict / std(error_measured_predict);
    qqplot(error_measured_predict_normalized);
    grid on
    title(['Q-Q plot error (predicted - measured) Anm(' allxticklabels{nn} ')']);
end
%}

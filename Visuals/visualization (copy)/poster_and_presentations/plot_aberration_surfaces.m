clear all
close all

params = set_parameters_xy_NAT;

groundtruth.global = [ ...
-30; %1 %defocus
0.4; %2 %defocus
0.0; %3 %defocus
0.0; %4 %defocus
0.0; %5 %astigmatism
0.0; %6 %astigmatism
0.0; %7 %astigmatism
0.0; %8 %astigmatism
0.0; %9 %astigmatism
0.0; %10 %coma
0.0; %11 %coma
0.0; %12 %coma
0.0; %13 %coma
0.0; %14 %spherical
]; %groundtruth gammas.

theta.global = [ ...
-10; %1 %defocus
-0.09; %2 %defocus
0.0; %3 %defocus
0.0; %4 %defocus
0.0; %5 %astigmatism
0.0; %6 %astigmatism
0.0; %7 %astigmatism
0.0; %8 %astigmatism
0.0; %9 %astigmatism
0.0; %10 %coma
0.0; %11 %coma
0.0; %12 %coma
0.0; %13 %coma
0.0; %14 %spherical
]; %fitted gammas.

nr_of_zernikes = size(params.aberrations,1);
orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernikes,1);
for jzer = 1:nr_of_zernikes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

xim = (params.imgSizeX*params.pixelsize*1e-3)*linspace(0,1,50);
yim = (params.imgSizeY*params.pixelsize*1e-3)*linspace(0,1,50);
[Xim,Yim] = meshgrid(xim,yim);

zernikeSurfaces_groundtruth = get_zernike_coefficients(Xim,Yim,groundtruth.global,params);
zernikeSurfaces_fitted = get_zernike_coefficients(Xim,Yim,theta.global,params);

zernikeSurfaces_groundtruth_normalized = 1e3*zernikeSurfaces_groundtruth/params.lambda;
zernikeSurfaces_fitted_normalized = 1e3*zernikeSurfaces_fitted/params.lambda;

zmin = min(cat(1,zernikeSurfaces_groundtruth_normalized,zernikeSurfaces_fitted_normalized),[],'all');
zmax = max(cat(1,zernikeSurfaces_groundtruth_normalized,zernikeSurfaces_fitted_normalized),[],'all');

% Plot surfaces per aberration
figure
    set(gcf,'Position',[25 90 1255 880])    
for nn = 1:nr_of_zernikes
    hsub = subplot(2,3,nn);    
    set(hsub, 'Visible', 'off');
    surf(Xim,Yim,zernikeSurfaces_groundtruth_normalized(:,:,nn),'FaceAlpha',0.4,'EdgeColor','none','FaceColor','#00F')
    hold on; axis square;
    surf(Xim,Yim,zernikeSurfaces_fitted_normalized(:,:,nn),'FaceAlpha',0.7,'EdgeColor','none','FaceColor','#F80')
    hold off; axis square;
    title(['A(' allxticklabels{nn} ')']);
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    zlabel('(m\lambda)');
    %zlim([zmin(nn) zmax(nn)]);
    zlim([zmin zmax]);
    ax=gca;
    ax.FontSize = 18;
end
%legend(hsub,'Groundtruth','Groundtruth - CRLB','Groundtruth + CRLB','Fitted surface','Location','east')
legend(hsub,'Groundtruth','Fitted surface','Location','east','Fontsize',16)
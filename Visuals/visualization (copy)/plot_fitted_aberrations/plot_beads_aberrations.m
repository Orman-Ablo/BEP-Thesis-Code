% Script to plot the aberrations found from z-stack beads over the FOV.
% A surface is interpolated between the beads.
% Author: Isabel Droste, TU Delft, 2024

% Input: theta, roixy and params

% Put here the number of pixels
params.imgSizeX = 512;
params.imgSizeY = 512;

%%%%%%%%%%%

pixelsize = params.pixelsize;

X = roixy(1,:)*pixelsize*1e-3;
Y = roixy(2,:)*pixelsize*1e-3;

xsize = pixelsize*params.imgSizeX;
ysize = pixelsize*params.imgSizeY;

zernike_coeffs = theta(6:end,:)*1e3/params.lambda; % zernike coefficients in mlambda
nr_of_zernikes = size(params.aberrations,1);
orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernikes,1);

for jzer = 1:nr_of_zernikes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

[Xq,Yq] = meshgrid(1:xsize*1e-3);

figure
for izer = 1:size(params.aberrations,1)
    hsub = subplot(3,5,izer);
    V = zernike_coeffs(izer,:);
    Vq = griddata(X,Y,V,Xq,Yq,"cubic");
    surf(Xq,Yq,Vq,'FaceAlpha',0.7)
    shading interp
    hold on
    plot3(X,Y,V,"o",'LineWidth',1.5,'MarkerSize',4)
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    zlabel('m\lambda')
    %zlim([-100 100]) % choose min and max of z-axis
    title(['A(' allxticklabels{izer} ')']);
    ax = gca;
    ax.FontSize = 14;

end
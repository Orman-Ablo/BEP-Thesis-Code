
%% Plot data.

% Plot individual spots
h = dipshow(allspots(:,:,:,:),'lin');
diptruesize(40000/params.Mx)
h.Name = 'allspots';
colormap parula

% Plot Zernike surfaces and spots (black dots)
%xim = (params.imgSizeX*params.pixelsize*1e-3)*linspace(0,1,50);
%yim = (params.imgSizeY*params.pixelsize*1e-3)*linspace(0,1,50);
xim = linspace(-1,1,50);
yim = linspace(-1,1,50);
[Xim,Yim] = meshgrid(xim,yim);

zernikeSurfaces = get_zernike_coefficients(Xim,Yim,groundtruth.global,params);

% Convert aberrations from nm to mlambda
zernikeCoefficients_normalized = 1e3*zernikeCoefficients/params.lambda;
zernikeSurfaces_normalized = 1e3*zernikeSurfaces/params.lambda;

zmin = min(zernikeSurfaces_normalized,[],"all");
zmax = max(zernikeSurfaces_normalized,[],"all");

if zmin==zmax
    zmin = zmin - 40;
    zmax = zmax + 40;
end

nr_of_zernikes = size(params.aberrations,1);

orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernikes,1);
for jzer = 1:nr_of_zernikes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

figure
for nn = 1:nr_of_zernikes
    subplot(2,3,nn)    
    surf(Xim,Yim,zernikeSurfaces_normalized(:,:,nn),'FaceAlpha',0.75)
    hold on; axis square; shading interp;
    scatter3(fov_coordinates(1,:),fov_coordinates(2,:),zernikeCoefficients_normalized(:,nn),3,'o','filled','MarkerFaceColor',[0 0 0]);
    title(['A(' allxticklabels{nn} ')']);
    xlabel('x (nm)');
    ylabel('y (nm)');
    %zlabel('Aberration coefficient (m\lambda)');
    zlabel('(m\lambda)');
    zlim([zmin zmax]);
    ax=gca;
    ax.FontSize = 16;
end

% Plot total wrms
wrms_surface = sqrt(sum(zernikeSurfaces_normalized.^2,3));
wrms_spots = sqrt(sum(zernikeCoefficients_normalized.^2,2));
figure
surf(Xim,Yim,wrms_surface,'FaceAlpha',0.75);
hold on; shading interp;
scatter3(fov_coordinates(1,:),fov_coordinates(2,:),wrms_spots,3,'o','filled','MarkerFaceColor',[0 0 0]);
title('Wrms');
xlabel('x (\mum)');
ylabel('y (\mum)');
zlabel('(m\lambda)');
ax=gca;
ax.FontSize = 18;

% Plot averaged PSFs
upsfac = 20;
imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;
Npatch = 4;

allpositions = groundtruth.local([2 1],:)/params.pixelsize;

figure
for xi = 1:Npatch
    for yi = 1:Npatch

        xrange = [floor((imgSizeX/Npatch)*(xi-1))+1,floor((imgSizeX/Npatch))*xi];
        yrange = [floor((imgSizeY/Npatch)*(yi-1))+1,floor((imgSizeY/Npatch))*yi];
        indices = find((roixy(1,:) >= xrange(1) & roixy(1,:) <= xrange(2) & roixy(2,:) >= yrange(1) & roixy(2,:) <= yrange(2)));
        indices = indices(indices<2e5);
        Ncfg = numel(indices);

        % select spots with these indices
        allspots_selected = allspots(:,:,:,indices);
        allpositions_selected = allpositions(:,indices);

        % feed spot stack to function for shift upsample and sum
        PSFsum = sum_shift_ups_PSFs(allspots_selected,allpositions_selected,upsfac);
        
        subplot(Npatch,Npatch,(xi-1)*Npatch+yi)
        imagesc(PSFsum)
        title(sprintf('(x=%i, y=%i) modeled',xi,yi))
        axis square
        axis off
        ax = gca;
        ax.FontSize = 8;

    end
end

% % plot single spot and emitter locations
% gt_x = groundtruth.local(1,:);
% gt_y = groundtruth.local(2,:);
% xrange = params.pixelsize*[-floor(params.Mx/2) floor(params.Mx/2)];
% yrange = params.pixelsize*[-floor(params.My/2) floor(params.My/2)];
% 
% figure
% imagesc(xrange,yrange, allspots(:,:,:,1))
% hold on
% scatter(gt_y(1),gt_x(1),100,'x','red')
% axis square;
% xlabel('y');
% ylabel('x');
% ax=gca;
% ax.FontSize = 20;

%     % Plot distribution of sum of pixel intensities
%     N_sums = squeeze(sum(allspots,1:3));
%     figure
%     histogram(N_sums,20)
%     ax=gca;
%     ax.FontSize = 20;
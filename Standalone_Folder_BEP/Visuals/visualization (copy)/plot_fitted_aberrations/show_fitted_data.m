%% Process

max_merit_index = 1;

theta.local = theta_full.local(:,:,max_merit_index);
theta.global = theta_full.global(:,:,max_merit_index);
thetastore.local = thetastore_full.local(:,:,:,max_merit_index);
thetastore.global = thetastore_full.global(:,:,max_merit_index);
meritstore.local = meritstore_full.local(:,:,max_merit_index);
meritstore.global = meritstore_full.global(:,:,max_merit_index);
mu = mu_full(:,:,:,:,max_merit_index);
dmudtheta.local = dmudtheta_full.local(:,:,:,:,:,max_merit_index);
dmudtheta.global = dmudtheta_full.global(:,:,:,:,:,max_merit_index);
outliers = outliers_full(max_merit_index);
outliers = outliers{1};
no_outliers = setdiff(1:params.Ncfg,outliers);

fitmodel = params.fitmodel;

% Mean and standard deviation of error
if params.groundtruth_exists
    error_local = theta.local - groundtruth.local;
    error_xy = error_local(1:2,no_outliers);
    
    mean_position_error_xy = mean(error_xy,2)';
    
    distance_from_mean = sqrt((error_xy(1,:)-mean_position_error_xy(1)).^2 + (error_xy(2,:)-mean_position_error_xy(2)).^2);
    std_error_xy = sqrt(mean(distance_from_mean.^2));
    
    if strcmp(fitmodel,'xyz-gamma')
        error_local_z = error_local(3,no_outliers);
        mean_error_z = mean(error_local(3,no_outliers));
        std_error_z = std(error_local(3,no_outliers));
        error_local_N = error_local(4,no_outliers);
        mean_error_N = mean(error_local(4,no_outliers));
        std_error_N = std(error_local(4,no_outliers));
        error_local_b = error_local(5,no_outliers);
        mean_error_b = mean(error_local(5,no_outliers));
        std_error_b = std(error_local(5,no_outliers));

    elseif strcmp(fitmodel,'xy-gamma')
        error_local_N = error_local(3,:);
        mean_error_N = mean(error_local(3,:));
        std_error_N = std(error_local(3,:));
        error_local_b = error_local(4,:);
        mean_error_b = mean(error_local(4,:));
        std_error_b = std(error_local(4,:));
    end
end

% Compute crlb local.
[crlb,rcondstore_local,Fisher_local] = get_fisher_crlb(params,mu,dmudtheta.local);

% Average CRLB
crlb_avg = mean(crlb,2);
crlb_avg_xy = sqrt(2)*mean([crlb_avg(1) crlb_avg(2)]);
if strcmp(fitmodel,'xyz-gamma')
    crlb_avg_z = crlb_avg(3);
    crlb_avg_N = crlb_avg(4);
    crlb_avg_b = crlb_avg(5);
elseif strcmp(fitmodel,'xy-gamma')
    crlb_avg_N = crlb_avg(3);
    crlb_avg_b = crlb_avg(4);
end

xn = linspace(-1,1,50);
yn = linspace(-1,1,50);
[Xn,Yn] = meshgrid(xn,yn);

%%
% Compute fit error
errM = get_fiterror(mu,allspots,params);

% Compute average Wrms value
if params.groundtruth_exists
    [avg_wrms_error,wrms_error] = get_wrms_error(theta,groundtruth,roixy,params,true);
end

%% Local parameter plots

% Position error x,y
if params.groundtruth_exists
    figure;
    points_xy = scatter(error_local(1,:),error_local(2,:),200,'Marker','.','MarkerEdgeColor','black');
    hold on
    mean_error = scatter(mean_position_error_xy(1),mean_position_error_xy(2),400,'Marker','x','MarkerEdgeColor','blue','LineWidth',4);
    hold on
    std_data = viscircles(mean_position_error_xy,std_error_xy,'Color','blue','LineWidth',4);
    hold on
    zero_error_xy = scatter(0,0,400,'Marker','x','MarkerEdgeColor','red','LineWidth',4);
    hold on
    bound_xy = viscircles([0 0],crlb_avg_xy,'Color','r','LineWidth',4);
    hold on

    axis equal
    legend([points_xy std_data bound_xy],{'error x,y','\sigma(error x,y)', 'CRLB'});
    xlabel('x (nm)','FontSize',22);
    ylabel('y (nm)','FontSize',22);
    title('Position error x,y');
    ax=gca;
    ax.FontSize = 20;

    % Position error z
    if strcmp(fitmodel,'xyz-gamma')
        ax = axes(figure, 'NextPlot', 'add', 'XColor', 'none');
        y_axis = error_local_z;
        x_axis = 0.5*crlb_avg_z*randn(numel(y_axis),1);
        points_z = scatter(x_axis,y_axis,75,'Marker','.','MarkerEdgeColor','black');
        std_z = errorbar(0,mean_error_z,std_error_z,'LineWidth',3,'MarkerSize',10,'Marker','x','MarkerEdgeColor','blue','Color','blue','CapSize',50);
        bound_z = errorbar(0,0,crlb_avg_z,'LineWidth',3,'MarkerSize',10,'Marker','x','MarkerEdgeColor','red','Color','red','CapSize',50);

        axis equal
        legend([points_z std_z bound_z],{'error z', '\sigma(error z)','CRLB'});
        ylabel('z (nm)','FontSize',20);
        title('Position error z');
        ax=gca;
        ax.FontSize = 20;
    end

    % Photon count error N
    ax = axes(figure, 'NextPlot', 'add', 'XColor', 'none');
    y_axis = error_local_N;
    x_axis = 0.5*crlb_avg_N*randn(numel(y_axis),1);
    points_N = scatter(x_axis,y_axis,75,'Marker','.','MarkerEdgeColor','black');
    std_N = errorbar(0,mean_error_N,std_error_N,'LineWidth',3,'MarkerSize',10,'Marker','x','MarkerEdgeColor','blue','Color','blue','CapSize',50);
    bound_N = errorbar(0,0,crlb_avg_N,'LineWidth',3,'MarkerSize',10,'Marker','x','MarkerEdgeColor','red','Color','red','CapSize',50);

    axis equal
    legend([points_N std_N bound_N],{'error N', '\sigma(error N)','CRLB'});
    ylabel('N(counts)','FontSize',20);
    title('Photon count error N');
    ax=gca;
    ax.FontSize = 20;

    % Background error b
    ax = axes(figure, 'NextPlot', 'add', 'XColor', 'none');
    y_axis = error_local_b;
    x_axis = 0.5*crlb_avg_b*randn(numel(y_axis),1);
    points_b = scatter(x_axis,y_axis,75,'Marker','.','MarkerEdgeColor','black');
    std_b = errorbar(0,mean_error_b,std_error_b,'LineWidth',3,'MarkerSize',10,'Marker','x','MarkerEdgeColor','blue','Color','blue','CapSize',50);
    bound_b = errorbar(0,0,crlb_avg_b,'LineWidth',3,'MarkerSize',10,'Marker','x','MarkerEdgeColor','red','Color','red','CapSize',50);

    axis equal
    legend([points_b std_b bound_b],{'error b', '\sigma(error b)','CRLB'});
    ylabel('b(counts per pixel)','FontSize',20);
    title('Background error b');
    ax=gca;
    ax.FontSize = 20;
end

%% Plot surfaces

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

[Xim,Yim] = normalized_to_physical_fov_coordinates(Xn,Yn,params);

% %% Plot total Wrms
% figure('Name','Total Wrms')
% if params.groundtruth_exists
%     surf(Xim,Yim,wrms_surface_groundtruth,'FaceAlpha',0.4,'EdgeColor','none','FaceColor','#00F');
%     hold on
% end
% surf(Xim,Yim,wrms_surface_fitted,'FaceAlpha',0.7,'EdgeColor','none','FaceColor','#F80');
% title('Wrms');
% xlabel('x (\mum)');
% ylabel('y (\mum)');
% zlabel('(m\lambda)');
% ax=gca;
% ax.FontSize = 18;
% if params.groundtruth_exists
%     legend('Groundtruth','Fitted surface','Location','east','Fontsize',22)
% else
%     legend('Fitted surface','Location','east','Fontsize',22)
% end

%% Plot aberration surfaces
figure('Name','Aberration surfaces')  
for nn = 1:nr_of_zernikes
    hsub = subplot(2,3,nn);    
    set(hsub, 'Visible', 'off');
    if params.groundtruth_exists
        surf(Xim,Yim,zernikeSurfaces_groundtruth_normalized(:,:,nn),'FaceAlpha',0.4,'EdgeColor','none','FaceColor','#00F')
        hold on; axis square;
    end
    surf(Xim,Yim,zernikeSurfaces_fitted_normalized(:,:,nn),'FaceAlpha',0.7,'EdgeColor','none','FaceColor','#F80')
    hold off; axis square;
    title(['A(' allxticklabels{nn} ')']);
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    zlabel('(m\lambda)');
    ax=gca;
    ax.FontSize = 18;
end
%legend(hsub,'Groundtruth','Groundtruth - CRLB','Groundtruth + CRLB','Fitted surface','Location','east')
if params.groundtruth_exists
    legend(hsub,'Groundtruth','Fitted surface','Location','east','Fontsize',22)
else
    legend(hsub,'Fitted surface','Location','east','Fontsize',22)
end


%% Plot error surfaces and CRLB
% 
% [CRLB_global_surface,Fisher_global,Fisher_global_inverse,CRLB_gammas,CRLB_Wsys,rcondstore_global] = get_fisher_crlb_global(mu,dmudtheta.global,Xn,Yn,outliers,theta.global,params);
% CRLB_global_surface = 1e3*CRLB_global_surface/params.lambda;
% 
% if params.groundtruth_exists
%     error_surfaces = abs(zernikeSurfaces_groundtruth_normalized - zernikeSurfaces_fitted_normalized);
%     zmin = min(cat(1,error_surfaces,CRLB_global_surface),[],[1 2]);
%     zmax = max(cat(1,error_surfaces,CRLB_global_surface),[],[1 2]);
% else 
%     zmin = min(CRLB_global_surface,[],[1 2]);
%     zmax = max(CRLB_global_surface,[],[1 2]);
% end
% 
% figure('Name','Aberration CRLB')
% for nn = 1:nr_of_zernikes
%     hsub = subplot(3,3,nn);    
%     set(hsub, 'Visible', 'off');
%     if params.groundtruth_exists
%         surf(Xim,Yim,error_surfaces(:,:,nn),'FaceAlpha',0.4,'EdgeColor','none','FaceColor','#00F')
%         hold on; axis square;
%     end
%     surf(Xim,Yim,CRLB_global_surface(:,:,nn),'FaceAlpha',0.7,'EdgeColor','none','FaceColor','#F80')
%     hold off; axis square;
%     title(['A(' allxticklabels{nn} ')']);
%     xlabel('x (\mum)');
%     ylabel('y (\mum)');
%     zlabel('(m\lambda)');
%     if zmin(nn)==zmax(nn)
%         zmin(nn) = zmin(nn) - 10;
%         zmax(nn) = zmax(nn) + 10;
%     end
%     zlim([zmin(nn) zmax(nn)]);
%     %zlim([zmin zmax]);
%     ax=gca;
%     ax.FontSize = 18;
% end
% 
% if params.groundtruth_exists
%     legend(hsub,'Error with groundtruth','CRLB','Location','east','Fontsize',16)
% else
%     legend(hsub,'CRLB','Location','east','Fontsize',16)
% end

% %% Plot Wrms error surface
% 
% if params.groundtruth_exists
%     wrms_error_surface = sqrt(sum((error_surfaces).^2,3));
%     crlb_surface_total = sqrt(sum((CRLB_global_surface).^2,3));
%     figure('Name','Wrms error')
%     surf(Xim,Yim,wrms_error_surface,'FaceAlpha',0.45,'EdgeColor','none','FaceColor','#00F');
%     hold on;
%     surf(Xim,Yim,crlb_surface_total,'FaceAlpha',0.6,'EdgeColor','none','FaceColor','#F80');
%     title('Wrms error');
%     xlabel('x (\mum)');
%     ylabel('y (\mum)');
%     zlabel('(m\lambda)');
%     ax=gca;
%     ax.FontSize = 18;
%     legend('Wrms error','CRLB','Location','east','Fontsize',16)
% end

%% Plot measured vs. fitted spots.
measured_fitted = cat(2,allspots,mu);
dipshow(measured_fitted(:,:,:,no_outliers),'lin','name','Measured ROIs vs. modeled PSFs (without outliers)');
diptruesize(40000/params.Mx)
colormap parula

% dipshow(measured_fitted(:,:,:,outliers),'lin','name','outliers');
% diptruesize(40000/params.Mx)
% colormap parula

%% Averaged and upsampled PSFs

upsfac = 10;
imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;
Npatch = 3;
zrange = [-1000,1000];
z_positions = theta.local(3,no_outliers);

allpositions = theta.local([2 1],:)/params.pixelsize;

figure
for xi = 1:Npatch
    for yi = 1:Npatch

        xrange = [floor((imgSizeX/Npatch)*(xi-1))+1,floor((imgSizeX/Npatch))*xi];
        yrange = [floor((imgSizeY/Npatch)*(yi-1))+1,floor((imgSizeY/Npatch))*yi];
        indices = find((roixy(1,no_outliers) >= xrange(1) & roixy(1,no_outliers) <= xrange(2) & roixy(2,no_outliers) >= yrange(1) & roixy(2,no_outliers) <= yrange(2) ...
            & z_positions > zrange(1) & z_positions < zrange(2) ));
        indices = indices(indices<2e5);
        Ncfg = numel(indices);

        % select spots with these indices
        allspots_selected = allspots(:,:,:,indices);
        mu_selected = mu(:,:,:,indices);
        allpositions_selected = allpositions(:,indices);

        % feed spot stack to function for shift upsample and sum
        PSFsum_measured = sum_shift_ups_PSFs(allspots_selected,allpositions_selected,upsfac);
        PSFsum_modeled= sum_shift_ups_PSFs(mu_selected,allpositions_selected,upsfac);
        
        subplot(Npatch,Npatch*2,(xi-1)*2*Npatch+2*yi-1);
        imagesc(PSFsum_measured)
        title(sprintf('Measured PSF (x=%i, y=%i) modeled',xi,yi))
        axis square
        axis off
        ax = gca;
        ax.FontSize = 8;

        subplot(Npatch,Npatch*2,(xi-1)*2*Npatch+2*yi);
        imagesc(PSFsum_modeled)
        title(sprintf('Modeled PSF (x=%i, y=%i) modeled',xi,yi))
        axis square
        axis off
        ax = gca;
        ax.FontSize = 8;

    end
end


%% plot gammas vs. iterations
% %iter_max = 40;
% gamma_labels = {'gamma1', 'gamma2', 'gamma3', 'gamma4', 'gamma5', 'gamma6', 'gamma7', 'gamma8', 'gamma9', 'gamma10', 'gamma11', 'gamma12', 'gamma13','gamma14','gamma15','gamma16'};
% iters = 1:size(thetastore.global,2);
% %iters = 1:iter_max;
% figure('name', 'gammas vs. iterations')
% for g=1:params.numgammas
%     subplot(3,6,g);
%     gamma_iter = thetastore.global(g,:);
%     gamma_iter_plot = gamma_iter;
%     %gamma_iter_plot = gamma_iter(1:iter_max);
%     if params.groundtruth_exists
%         plot(iters,groundtruth.global(g)*ones(size(iters)),'LineWidth',4);
%         hold on
%     end
%     plot(iters,gamma_iter_plot,'LineWidth',3);
%     xlabel('iterations');
%     ylabel(gamma_labels{g});
%     ax=gca;
%     ax.FontSize = 14;
%     ylim("padded");
% end
% 
% if params.groundtruth_exists
%     legend('Groundtruth','gamma');
% else
%     legend('gamma');
% end

%% Scatter plot of emitters
% figure;
% scatter3(localizations(:,2),localizations(:,3),localizations(:,4))
% xlabel('x (nm)')
% ylabel('y (nm)')
% zlabel('z (nm)')

% Plot average and std aberration surfaces and compare to CRLB

close all
clear all

path_input_data = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Yutong_Nanorulers\data_analysis\nanorulers\fitted_aberrations\fitted_aberrations_ROI7_NR80R_1114_Oil100x1-5_Col200_647P200-100_80ms_EM100_10000f_FL001_stk_1_*';

params = set_parameters_yutong_nanorulers;

nr_of_zernike_modes = size(params.aberrations,1);

all_input_data = dir(path_input_data);
nr_of_input_data = length(all_input_data);
fprintf('\n Found %i matching output data files\n', nr_of_input_data);

grid_size = 50;
xn = linspace(-1,1,grid_size);
yn = linspace(-1,1,grid_size);
[Xn,Yn] = meshgrid(xn,yn);

aberration_surfaces = zeros(nr_of_input_data,grid_size,grid_size,nr_of_zernike_modes);
aberration_CRLB_surfaces = zeros(nr_of_input_data,grid_size,grid_size,nr_of_zernike_modes);

%% Load data
for i = 1:nr_of_input_data
    disp(all_input_data(i).name);
    input_path = fullfile(all_input_data(i).folder,all_input_data(i).name);
    load(input_path);

    theta_global = theta_full.global(:,:,max_merit_index);
    aberration_surfaces(i,:,:,:) = 1e3*get_zernike_coefficients(Xn,Yn,theta_global,params)/params.lambda;

    mu = mu_full(:,:,:,:,max_merit_index);
    dmudtheta.global = dmudtheta_full.global(:,:,:,:,:,max_merit_index);
    aberration_CRLB_surfaces(i,:,:,:) = 1e3*get_fisher_crlb_global(mu,dmudtheta.global,Xn,Yn,params)/params.lambda;

end

aberration_surfaces_mean = squeeze(mean(aberration_surfaces,1));
aberration_surfaces_std = squeeze(std(aberration_surfaces,1));
aberration_CRLB_surfaces_mean = squeeze(mean(aberration_CRLB_surfaces,1));

%% Plot mean surfaces
[Xim,Yim] = normalized_to_physical_fov_coordinates(Xn,Yn,params);

orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernike_modes,1);
for jzer = 1:nr_of_zernike_modes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

zmin = min(aberration_surfaces_mean,[],'all');
zmax = max(aberration_surfaces_mean,[],'all');

figure
for nn = 1:nr_of_zernike_modes
    hsub = subplot(2,3,nn);    
    set(hsub, 'Visible', 'off');
    surf(Xim,Yim,aberration_surfaces_mean(:,:,nn),'FaceAlpha',0.6,'EdgeColor','none','FaceColor','#F80')
    hold on; axis square;
    title(['A(' allxticklabels{nn} ')']);
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    zlabel('(m\lambda)');
    zlim([zmin zmax]);
    ax=gca;
    ax.FontSize = 16;
end
legend('Mean aberration surface','Location','east','Fontsize',16)

%% Plot std surfaces
[Xim,Yim] = normalized_to_physical_fov_coordinates(Xn,Yn,params);

orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernike_modes,1);
for jzer = 1:nr_of_zernike_modes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

zmin = min(cat(3,aberration_surfaces_std,aberration_CRLB_surfaces_mean),[],'all');
zmax = max(cat(3,aberration_surfaces_std,aberration_CRLB_surfaces_mean),[],'all');

figure
for nn = 1:nr_of_zernike_modes
    hsub = subplot(2,3,nn);    
    set(hsub, 'Visible', 'off');
    surf(Xim,Yim,aberration_surfaces_std(:,:,nn),'FaceAlpha',0.6,'EdgeColor','none','FaceColor','#F80')
    hold on; axis square;
    surf(Xim,Yim,aberration_CRLB_surfaces_mean(:,:,nn),'FaceAlpha',0.4,'EdgeColor','none','FaceColor','#00F')
    title(['A(' allxticklabels{nn} ')']);
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    zlabel('(m\lambda)');
    zlim([zmin zmax]);
    ax=gca;
    ax.FontSize = 16;
end
legend('Standard deviation','CRLB','Location','east','Fontsize',16)

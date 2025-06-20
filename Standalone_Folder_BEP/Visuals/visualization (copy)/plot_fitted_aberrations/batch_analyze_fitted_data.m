close all
clear all

%path_input_data = 'C:\Users\idroste\Desktop\TUDelft\Data\NAT_data\sample_data\input_spots_1.mat';
%path_input_data = 'C:\Users\idroste\Desktop\TUDelft\Data\NAT_data\experiments_MBL2023\simulation_20230324_063850\input_spots_*.mat';
%path_input_data = 'U:\SMLMaberrations\NATaberrations\SimulationData\simulation_20230418_070529\input_spots_*.mat';
%path_input_data = 'N:\tnw\IST\QI\users\idroste\Data\SimulationData\simulation_20230418_070529\input_spots_*.mat';


%path_output_data = "C:\Users\idroste\Desktop\TUDelft\Data\NAT_data\results\results_1.mat";
%path_output_data = 'C:\Users\idroste\Desktop\TUDelft\Data\NAT_data\experiments_MBL2023\simulation_20230324_063850\results_*.mat';
%path_output_data = 'U:\SMLMaberrations\NATaberrations\SimulationData\simulation_20230418_070529\results_*.mat';
%path_output_data = 'N:\tnw\IST\QI\users\idroste\Data\SimulationData\simulation_20230418_070529\results_*.mat';
%path_output_data = 'U:\isabeldrostePhD\ExperimentalData\DatasetShareLoc_Singh\results\NvsVariance\results_*.mat';
%path_output_data = 'U:\isabeldrostePhD\ExperimentalData\DatasetShareLoc_Singh\results\NvsVariance\500\results_*.mat';
%path_output_data = 'U:\isabeldrostePhD\ExperimentalData\DatasetShareLoc_Singh\results\results_linuxpc\results_*.mat';
path_output_data = 'C:\Users\idroste\Desktop\TUDelft\Data\NAT_data\LowPhotonCounts\2500\results\results_*';

%all_input_data = dir(path_input_data);
%nr_of_runs = length(all_input_data);
%fprintf('\n Found %i matching input data files\n', length(all_input_data));
% for i = 1:nr_of_runs
%     disp(all_input_data(i).name);
% end

all_output_data = dir(path_output_data);
nr_of_output_data = length(all_output_data);
nr_of_output_data = 10;
fprintf('\n Found %i matching output data files\n', length(all_output_data));
for i = 1:length(all_output_data)
    disp(all_output_data(i).name);
end

% aberration_threshold = 20;
% maximum_was_found = zeros(nr_of_runs,1);
% nr_of_maxima_found = zeros(nr_of_runs,1);
% nr_of_flipped_maxima_found = zeros(nr_of_runs,1);
% crlb_all = 0;
% location_bias_aberrations_all = [];
% location_bias_no_aberrations_all = [];
% location_std_aberrations_all = [];
% location_std_no_aberrations_all = [];
% wrms_all = [];
% spherical_max = [];
% spherical_nomax = [];
% indices_nomax = [];

params = set_parameters_xy_NAT;
nr_of_zernike_modes = size(params.aberrations,1);


xn = linspace(-1,1,50);
yn = linspace(-1,1,50);
[Xn,Yn] = meshgrid(xn,yn);

size_grid = size(Xn);
aberration_surfaces = zeros(nr_of_output_data,size_grid(1),size_grid(2),nr_of_zernike_modes);

for r=1:nr_of_output_data
    r
    output_path = fullfile(all_output_data(r).folder,all_output_data(r).name);
    load(output_path);
    
    %theta_global = theta_full.global(:,:,1).*[-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,-1]';
    theta_global = theta_full.global(:,:,1);
    aberration_surfaces(r,:,:,:) = 1e3*get_zernike_coefficients(Xn,Yn,theta_global,params)/params.lambda;
end


%% Plot

[Xim,Yim] = normalized_to_physical_fov_coordinates(Xn,Yn,params);

std_surfaces = squeeze(std(aberration_surfaces));

orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernike_modes,1);
for jzer = 1:nr_of_zernike_modes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

figure
for nn = 1:nr_of_zernike_modes
    hsub = subplot(2,3,nn);    
    set(hsub, 'Visible', 'off');
    surf(Xim,Yim,std_surfaces(:,:,nn),'FaceAlpha',0.4,'EdgeColor','none','FaceColor','#00F')
    hold on; axis square;
    title(['A(' allxticklabels{nn} ')']);
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    zlabel('(m\lambda)');
    %zlim([zmin(nn) zmax(nn)]);
    %zlim([zmin zmax]);
    ax=gca;
    ax.FontSize = 18;
end
sgtitle('Standard deviation surfaces')

%%

output_path = fullfile(all_output_data(1).folder,all_output_data(1).name);
load(output_path);

max_merit_index = 1;
theta.local = theta_full.local(:,:,max_merit_index);
theta.global = theta_full.global(:,:,max_merit_index);
%theta.global = theta.global.*[-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,-1]';
thetastore.local = thetastore_full.local(:,:,:,max_merit_index);
thetastore.global = thetastore_full.global(:,:,max_merit_index);
%meritstore = meritstore_full(:,:,max_merit_index);
%alambdastore.local = alambdastore_full.local(:,:,max_merit_index);
%alambdastore.global = alambdastore_full.global(:,:,max_merit_index);
%meritstore_total = meritstore_total_full(:,:,max_merit_index);
mu = mu_full(:,:,:,:,max_merit_index);
dmudtheta.local = dmudtheta_full.local(:,:,:,:,:,max_merit_index);
dmudtheta.global = dmudtheta_full.global(:,:,:,:,:,max_merit_index);

[CRLB_global_surface,Fisher_global,Fisher_global_inverse,rcondstore_global,] = get_fisher_crlb_global(mu,dmudtheta.global,Xn,Yn,params);
CRLB_global_surface = 1e3*CRLB_global_surface/params.lambda;

zmin_separate = min(cat(1,std_surfaces,CRLB_global_surface),[],[1 2]);
zmax_separate = max(cat(1,std_surfaces,CRLB_global_surface),[],[1 2]);

zmin = min(cat(1,std_surfaces(:,:,2:end),CRLB_global_surface(:,:,2:end)),[],'all');
zmax = max(cat(1,std_surfaces(:,:,2:end),CRLB_global_surface(:,:,2:end)),[],'all');

figure   
for nn = 1:nr_of_zernike_modes
    hsub = subplot(2,3,nn);    
    set(hsub, 'Visible', 'off');
    surf(Xim,Yim,CRLB_global_surface(:,:,nn),'FaceAlpha',0.6,'EdgeColor','none','FaceColor','#F80')
    hold on; axis square;
    surf(Xim,Yim,std_surfaces(:,:,nn),'FaceAlpha',0.4,'EdgeColor','none','FaceColor','#00F')
    title(['A(' allxticklabels{nn} ')']);
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    zlabel('(m\lambda)');
    zlim([zmin_separate(nn) zmax_separate(nn)]);
    ax=gca;
    ax.FontSize = 18;
    ztickformat('%.1f');
end
legend('CRLB surfaces','Standard deviation','Location','east','Fontsize',16)

%% Calculate RMS CRLB and variance

crlb_surface_total = sqrt(sum((CRLB_global_surface).^2,3));
wrms_std_surface = sqrt(sum((std_surfaces).^2,3));

average_crlb = mean(crlb_surface_total,"all")
average_std = mean(wrms_std_surface,"all")

max_crlb = max(crlb_surface_total,[],"all")
max_std = max(wrms_std_surface,[],"all")

%%
% for r=1:nr_of_output_data
%     r
%     apply_flip = false;
% 
%     output_path = fullfile(all_output_data(r).folder,all_output_data(r).name);
%     load(output_path);
%     
%     split_name = split(all_output_data(r).name,'_');
%     rr = str2num(split_name{2});
% 
%     if params.fit_aberrations
% 
%         % compute CRLB locations
%         
%         
%         nr_of_maxima_found(rr) = sum(wrms_error_average(:,1)<aberration_threshold);
%         nr_of_flipped_maxima_found(rr) = sum(wrms_error_average(:,2)<aberration_threshold);
%         if min(wrms_error_average(max_merit_index,:))<aberration_threshold
%             maximum_was_found(rr) = 1;
%         end
%             % other measure for is maximum was found.
% %         if nr_of_maxima_found(rr)
% %             maximum_was_found(rr) = 1;
% %         end
% 
%         if wrms_error_average(max_merit_index,1)>wrms_error_average(max_merit_index,2)
%             apply_flip = true;
%         end
%         
%         % bias
%         if maximum_was_found(rr)
%             mu = mu_(:,:,:,:,max_merit_index);
%             dmudtheta_local = dmudtheta_.local(:,:,:,:,:,max_merit_index);
%             crlb_avg = mean(get_fisher_crlb(params,mu,dmudtheta_local,'local'),2);
%             crlb_avg_xy = sqrt(2)*mean([crlb_avg(1) crlb_avg(2)]);
%             crlb_all = crlb_all + crlb_avg_xy;
%             
%             location_bias_aberrations_all(end+1) = location_bias(max_merit_index);
%             location_std_aberrations_all(end+1) = location_std(max_merit_index);
% 
%             wrms_all(end+1) = wrms_error_average(max_merit_index,apply_flip+1);
%             spherical_max(end+1) = groundtruth.global(14);
%         else
%             spherical_nomax(end+1) = groundtruth.global(14);
%             indices_nomax(end+1) = rr;
%         end
% 
%     else
% 
%         % bias
%         if maximum_was_found(rr)
%             location_bias_no_aberrations_all(end+1) = location_bias(max_merit_index);
%             location_std_no_aberrations_all(end+1) = location_std(max_merit_index);
%         end
% 
%     end
% 
% 
% end
% %%
% total_nr_of_maxima_found = sum(maximum_was_found);
% average_bias_aberrations = mean(location_bias_aberrations_all);
% average_bias_no_aberrations = mean(location_bias_no_aberrations_all);
% average_variance_aberrations = mean(location_std_aberrations_all)/(crlb_all/total_nr_of_maxima_found);
% average_variance_no_aberrations = mean(location_std_no_aberrations_all)/(crlb_all/total_nr_of_maxima_found);
% wrms_average = mean(wrms_all);
% 
% 
% %% Plot histograms
% 
% data = wrms_all;
% histogram(data);
% xline(mean(data), 'Color', 'r', 'LineWidth', 4);
% xlabel('Wrms error  (m\lambda)');
% ylabel('Counts');
% title('Wrms error');
% ax=gca;
% ax.FontSize = 28;

% data = location_bias_aberrations_all;
% histogram(location_bias_aberrations_all,10)
% xline(mean(data), 'Color', 'r', 'LineWidth', 4);
% xlabel('Location bias (nm)');
% ylabel('Counts');
% ax=gca;
% ax.FontSize = 24;
% 
% data = location_bias_no_aberrations_all;
% histogram(location_bias_no_aberrations_all,10)
% xline(mean(data), 'Color', 'r', 'LineWidth', 4);
% xlabel('Location bias (nm)');
% ylabel('Counts');
% ax=gca;
% ax.FontSize = 24;
% 
% data = location_std_aberrations_all;
% histogram(location_std_aberrations_all,8)
% xline(mean(data), 'Color', 'r', 'LineWidth', 4);
% xlabel('Standard deviation (nm)');
% ylabel('Counts');
% ax=gca;
% ax.FontSize = 24;
% 
% data = location_std_no_aberrations_all;
% histogram(location_std_no_aberrations_all,8)
% xline(mean(data), 'Color', 'r', 'LineWidth', 4);
% xlabel('Standard deviation (nm)');
% ylabel('Counts');
% ax=gca;
% ax.FontSize = 24;
% 
% data = 1e3*abs(spherical_nomax)/params.lambda;
% histogram(data,10)
% %xline(mean(data), 'Color', 'r', 'LineWidth', 2);
% xlabel('Spherical aberration (m\lambda)');
% ylabel('Counts');
% ax=gca;
% ax.FontSize = 24;


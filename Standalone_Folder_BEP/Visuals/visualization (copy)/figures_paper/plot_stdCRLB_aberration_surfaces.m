close all
clear all

path_output_data = '/home/idroste/Desktop/TUDelft/Data/Fu_fig4/data_analysis/fitted_aberrations/calculate_precision/005000spots/fitted_aberrations_*';

params = set_parameters_Fu_3D;

all_output_data = dir(path_output_data);
nr_of_output_data = length(all_output_data);
fprintf('\n Found %i matching data files\n', length(all_output_data));
for i = 1:nr_of_output_data
    disp(all_output_data(i).name);
end 

xn = linspace(-1,1,50);
yn = linspace(-1,1,50);
[Xn,Yn] = meshgrid(xn,yn);

nr_of_zernikes = size(params.aberrations,1);
orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernikes,1);

for jzer = 1:nr_of_zernikes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

zernike_surfaces = zeros(size(xn,2),size(yn,2),nr_of_zernikes,nr_of_output_data);

for r=1:nr_of_output_data
    r
    output_path = fullfile(all_output_data(r).folder,all_output_data(r).name);
    load(output_path);

    gammas_fitted = theta_full.global(:,:,max_merit_index);
    zernike_surfaces(:,:,:,r) = 1e3*get_zernike_coefficients(Xn,Yn,gammas_fitted,params)/params.lambda;
end

zmin = min(zernike_surfaces(:,:,[1 2 4 5 6],:),[],'all');
zmax = max(zernike_surfaces(:,:,[1 2 4 5 6],:),[],'all');
zmin_astigmatism = min(min(zernike_surfaces(:,:,3),[],'all'),0);
zmax_astigmatism = max(zernike_surfaces(:,:,3),[],'all');

%% Plot aberration surfaces
figure('Name','Aberration surfaces')  
for nn = 1:nr_of_zernikes
    hsub = subplot(2,3,nn);    
    set(hsub, 'Visible', 'off');
    surf(Xn,Yn,zernike_surfaces(:,:,nn,10),'FaceAlpha',0.7,'EdgeColor','none','FaceColor','#F80')
    hold off; axis square;
    title(['A(' allxticklabels{nn} ')']);
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    zlabel('(m\lambda)');
    %zlim([zmin(nn) zmax(nn)]);
    if nn == 3
        zlim([zmin_astigmatism zmax_astigmatism]);
    else
        zlim([zmin zmax]);
    end
    ax=gca;
    ax.FontSize = 18;
end
legend(hsub,'Fitted surface','Location','east','Fontsize',16)

%% Compute CRLB

r=1;
output_path = fullfile(all_output_data(r).folder,all_output_data(r).name);
load(output_path);

mu = mu_full(:,:,:,:,1);
dmudtheta.global = dmudtheta_full.global(:,:,:,:,:,1);
[CRLB_global_surfaces,Fisher_global,Fisher_global_inverse,CRLB_gammas,rcondstore_global] = get_fisher_crlb_global(mu,dmudtheta.global,Xn,Yn,outliers_full{1},params);
CRLB_global_surfaces_normalized = 1e3*CRLB_global_surfaces/params.lambda;


%% Plot std surfaces and CRLB

std_surfaces = std(zernike_surfaces,1,4);
avg_surfaces = mean(zernike_surfaces,4);

zmin_std = min(std_surfaces(:,:,[1 2 4 5 6],:),[],'all');
zmax_std = max(std_surfaces(:,:,[1 2 4 5 6],:),[],'all');
zmin_astigmatism_std = min(min(std_surfaces(:,:,3),[],'all'),0);
zmax_astigmatism_std = max(std_surfaces(:,:,3),[],'all');

figure('Name','Aberration surfaces')  
for nn = 1:nr_of_zernikes
    hsub = subplot(2,3,nn);    
    set(hsub, 'Visible', 'off');
    surf(Xn,Yn,std_surfaces(:,:,nn),'FaceAlpha',0.5,'EdgeColor','none','FaceColor','#00F')
    hold on; axis square;
    surf(Xn,Yn,CRLB_global_surfaces_normalized(:,:,nn),'FaceAlpha',0.5,'EdgeColor','none','FaceColor','#F80')
    title(['A(' allxticklabels{nn} ')']);
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    zlabel('(m\lambda)');
    %zlim([0 2.5]);
    % if nn == 3
    %     zlim([zmin_astigmatism_std zmax_astigmatism_std]);
    % else
    %     zlim([zmin_std zmax_std]);
    % end
    ax=gca;
    ax.FontSize = 18;
end
legend(hsub,'std surfaces','CRLB','Location','east','Fontsize',16)

%% Plot average surfaces and aber maps

load('U:\NATaberrations\ExperimentalData\DataNatureDeeplearningAberrations\demo2_FD_astig_NPC\demo2_FD_astig_NPC\aber_map.mat')
%%
aber_map_scaled = 1e3*aber_map;
% downsample
aber_map_downsampled = zeros(50, 50, 23);
for i = 1:23
    originalSlice = aber_map_scaled(:, :, i);
    downsampledSlice = imresize(originalSlice, [50, 50]);
    aber_map_downsampled(:, :, i) = downsampledSlice;
end

% rotate aber map
aber_map_adjusted = zeros(size(avg_surfaces));
aber_map_adjusted(:,:,1) = aber_map_downsampled(:,:,1);
aber_map_adjusted(:,:,3) = aber_map_downsampled(:,:,2);
aber_map_adjusted(:,:,4) = aber_map_downsampled(:,:,3);
aber_map_adjusted(:,:,5) = aber_map_downsampled(:,:,4);
aber_map_adjusted(:,:,6) = aber_map_downsampled(:,:,5);

aber_map_adjusted = permute(aber_map_adjusted,[2 1 3]);

zmin = min(cat(1,aber_map_adjusted,avg_surfaces),[],[1 2]);
zmax = max(cat(1,aber_map_adjusted,avg_surfaces),[],[1 2]);

figure
for nn=1:6
    hsub = subplot(2,3,nn);
    set(hsub, 'Visible', 'off');
    surf(Xn,Yn,aber_map_adjusted(:,:,nn),'FaceAlpha',0.5,'EdgeColor','none','FaceColor','#00F')
    hold on; axis square;
    surf(Xn,Yn,avg_surfaces(:,:,nn),'FaceAlpha',0.6,'EdgeColor','none','FaceColor','#F80')
    if zmin(nn)==zmax(nn)
        zmin(nn) = zmin(nn) - 1;
        zmax(nn) = zmax(nn) + 1;
    end
    zlim([zmin(nn) zmax(nn)]);
    %title(['A(' allxticklabels{nn} ')']);
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    zlabel('(m\lambda)');
    ax=gca;
    ax.FontSize = 18;
end

legend(hsub,'Aberration map paper','Average fitted surfaces','Location','east','Fontsize',16)

%% Plot standard deviation versus number of spots

clear all
path_output_data = '/home/idroste/Desktop/TUDelft/Data/Fu_fig4/data_analysis/fitted_aberrations/calculate_precision/*spots';

params = set_parameters_Fu_3D;

alldirs = dir(path_output_data);
ndirs = length(alldirs);
fprintf('\n Found %i directories\n', length(alldirs));

nr_of_zernikes = size(params.aberrations,1);
orders = params.aberrations(:,1:2);
allxticklabels = cell(nr_of_zernikes,1);
for jzer = 1:nr_of_zernikes
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end
    
xn = linspace(-1,1,50);
yn = linspace(-1,1,50);
[Xn,Yn] = meshgrid(xn,yn);
std_surfaces = zeros(50,50,nr_of_zernikes,ndirs);
CRLB_global_surfaces_dir = zeros(50,50,nr_of_zernikes,ndirs);
mean_surfaces = zeros(50,50,nr_of_zernikes,ndirs);

for d=1:23%ndirs
   
    path_output_data_dir = fullfile(alldirs(d).folder,alldirs(d).name)
    dirname = alldirs(d).name;

    all_output_data_dir = dir(path_output_data_dir);
    nr_of_output_data_dir = length(all_output_data_dir);
    fprintf('\n Found %i matching data files\n', length(all_output_data_dir));
    % for i = 1:nr_of_output_data_dir
    %     disp(all_output_data_dir(i).name);
    % end
    
    zernike_surfaces_dir = zeros(size(xn,2),size(yn,2),nr_of_zernikes,nr_of_output_data_dir-2);
    
    for r=3:nr_of_output_data_dir
        if(~all_output_data_dir(r).isdir)
            output_path = fullfile(all_output_data_dir(r).folder,all_output_data_dir(r).name);
    
            load(output_path);
        
            gammas_fitted = theta_full.global(:,:,max_merit_index);
            zernike_surfaces_dir(:,:,:,r-2) = 1e3*get_zernike_coefficients(Xn,Yn,gammas_fitted,params)/params.lambda;


            % CRLB
            mu = mu_full(:,:,:,:,1);
            dmudtheta.global = dmudtheta_full.global(:,:,:,:,:,1);
            [CRLB_global_surfaces,Fisher_global,Fisher_global_inverse,CRLB_gammas,rcondstore_global] = get_fisher_crlb_global(mu,dmudtheta.global,Xn,Yn,outliers_full{1},params);
            CRLB_global_surfaces_dir(:,:,:,d) = CRLB_global_surfaces_dir(:,:,:,d) + 1e3*CRLB_global_surfaces/params.lambda;
            
        end 
    end
    std_surfaces(:,:,:,d) = std(zernike_surfaces_dir,1,4);
    mean_surfaces(:,:,:,d) = mean(zernike_surfaces_dir,4);
    CRLB_global_surfaces_dir(:,:,:,d) = CRLB_global_surfaces_dir(:,:,:,d)/(nr_of_output_data_dir-2);
    
end

mean_std_surfaces = squeeze(mean(std_surfaces,[1 2]));
sum_mean_std_surfaces = sum(mean_std_surfaces,1);

mean_CRLB_surfaces = squeeze(mean(CRLB_global_surfaces_dir,[1 2]));
sum_CRLB_surfaces = sum(mean_CRLB_surfaces,1);

%%
mean_mean_surfaces = squeeze(mean(mean_surfaces,[1 2]));

%%
nspots = [75 100 200 300 400 500 600 700 800 900 1000 1250 1500 1750 2000 2500 3000 3500 4000 4500 5000 5500 6000];% 7000 8000 9000 10000];

figure
plot(nspots,sum_mean_std_surfaces(1:23),'-o','LineWidth',5,'MarkerSize',18,'Color',"#0072BD")
hold on
plot(nspots,sum_CRLB_surfaces(1:23),'-o','LineWidth',5,'MarkerSize',18,'Color',"#A2142F")
% hold on
% plot(nspots,350./sqrt(nspots),'-','LineWidth',2,'MarkerSize',7)
% hold on
% plot(nspots,70./sqrt(nspots),'-','LineWidth',2,'MarkerSize',7)
xlabel('M_s')
ylabel('Standard deviation (m\lambda)')
title('Sum of aberrations');
%ylim([0,12])
%legend('Standard deviation','CRLB')
ax=gca;
ax.FontSize = 40;

%%
nn=6;
figure
plot(nspots,mean_std_surfaces(nn,1:23),'-o','LineWidth',5,'MarkerSize',18,'Color',"#0072BD")
hold on
plot(nspots,mean_CRLB_surfaces(nn,1:23),'-o','LineWidth',5,'MarkerSize',18,'Color',"#A2142F")
xlabel('M_s')
ylabel('Standard deviation (m\lambda)')
title(['A(' allxticklabels{nn} ')']);
%ylim([0,12])
%legend('Standard deviation','CRLB')
ax=gca;
ax.FontSize = 40;
%%
nspots = [75 100 200 300 400 500 600 700 800 900 1000 1250 1500 1750 2000 2500 3000 3500 4000 4500 5000 5500 6000];% 7000 8000 9000 10000];

nn=4;
figure
plot(nspots,mean_mean_surfaces(nn,1:23),'-o','LineWidth',5,'MarkerSize',18,'Color',"#0072BD")
%hold on
%plot(nspots,mean_CRLB_surfaces(nn,1:23),'-o','LineWidth',5,'MarkerSize',18,'Color',"#A2142F")
xlabel('M_s')
ylabel('Standard deviation (m\lambda)')
title(['A(' allxticklabels{nn} ')']);
%ylim([0,12])
%legend('Standard deviation','CRLB')
ax=gca;
ax.FontSize = 40;

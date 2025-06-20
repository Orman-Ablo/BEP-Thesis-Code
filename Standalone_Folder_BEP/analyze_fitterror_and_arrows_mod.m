% 
% %mu = mu_beadaber;
% mu = mu_smlmaber;
% errM = get_fiterror(mu,allspots,params);
% mean(errM(1,:),2)
% std(errM(1,:),1,2)
% 
% %mean_chi_square = 1;
% mean_chi_square = params.Mx*params.My
% % std_chi_square = mean(sqrt(squeeze(2/(params.Mx*params.My)+mean(1./mu,1:3)./(params.Mx*params.My))))
% std_chi_square = mean(sqrt(squeeze(2*params.Mx*params.My+mean(1./mu,1:3))))
% 
% figure
% a = histogram(errM(1,:));
% binedges = a.BinEdges;
% bincounts = a.BinCounts;
% 
% %%
% errM_bead_aber = get_fiterror(mu_adapted_all_beads,allspots,params);
% 
% figure
% b = histogram(errM_bead_aber(1,:),BinEdges=binedges);
% bincounts_beads = b.BinCounts;
% 
% %%
% figure
% plot(binedges(1:end-1),bincounts,'LineWidth',2)
% hold on
% plot(binedges(1:end-1),bincounts_beads,'LineWidth',2)
% legend('from single molecules','from beads')

%%
%path_beadaber = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\Data_wo_Astigmatism\localization\localization_doubles\zero_aber\localization_mu_allspots\localization_mu_noaber_segmentation_roi9_thr20_Data0001.mat';
%path_smlmaber = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\Data_wo_Astigmatism\localization\localization_doubles\smlm_aber\localization_mu_allspots\localization_mu_segmentation_roi9_thr20_Data0001.mat';

path_beadaber = "C:\Users\Naam\Documents\BEP\data\final\localization_OTF_localization_segmentation_roi17_thr15p5_Data0001.mat";
path_smlmaber = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Extra input data\transfer_3055388_files_4e4cf1e2\localization_segmentation_roi17_thr15p5_Data0001.mat";

%path_beadaber = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\Data_w_Astigmatism\localization\localization_doubles\localizations_mu_allspots\localization_const_aber_mu_allspots\localization_segmentation_roi17_thr15p5_Data0001.mat';
%path_smlmaber = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\Data_w_Astigmatism\localization\localization_doubles\localizations_mu_allspots\localization_smlm_aber_mu_allspots\localization_segmentation_roi17_thr15p5_Data0001.mat';

data_second_fit = load(path_beadaber);
data_initial_fit = load(path_smlmaber);
%%


%% the code below is made to combine outlier position,s and isolate all but the correct locations. Maybe we can preemptively filter them?
%% filter the initial fit spots by using ID to find the indexes (works with combined files)
% (or by applying the outlier indexes(issue: they are file-dependent)
%% strategy: remove outliers from dataset, then sort by indexes of correct ID
params = data_second_fit.params;
roixy = data_second_fit.roixy_fit;
allspots = data_second_fit.Spots_fit;
mu_beadaber = data_initial_fit.mu;
mu_smlmaber = data_second_fit.mu;
localizations_2 = data_second_fit.localizations;
localizations_1 = data_initial_fit.localizations;


%%
outliers1 = data_initial_fit.outliers;
outliers2 = data_second_fit.outliers;

Ncfg2 = size(mu_smlmaber,4);
Ncfg1 = size(mu_beadaber,4);
no_outliers1 = setdiff(1:Ncfg1,outliers1);
no_outliers2 = setdiff(1:Ncfg2,outliers2);




Ncfg = size(roixy,2);
% no_outliers = setdiff(1:Ncfg,outliers_combined);

roixy = roixy(:,no_outliers2);
allspots = allspots(:,:,:,no_outliers2);

mu_smlmaber = mu_smlmaber(:,:,:,no_outliers2);



mu_beadaber = mu_beadaber(:,:,:,no_outliers1);

% in our case, the second file is a filter of the first file: we need to
% select the right values
[C,index_1,index_2] = intersect(localizations_1(:,1),localizations_2(:,1));
mu_beadaber = mu_beadaber(:,:,:,index_1);
z_positions = data_second_fit.localizations(:,4)';




%%
% Ncfg = size(roixy,2);
% outliers_beadaber = data_initial_fit.outliers;
% outliers_smlmaber = data_second_fit.outliers;
% outliers_combined = unique(cat(1,outliers_beadaber,outliers_smlmaber));
% no_outliers = setdiff(1:Ncfg,outliers_combined);
% 
% roixy = roixy(:,no_outliers);
% allspots = allspots(:,:,:,no_outliers);
% mu_beadaber = mu_beadaber(:,:,:,no_outliers);
% mu_smlmaber = mu_smlmaber(:,:,:,no_outliers);
% z_positions = data_second_fit.localizations_with_outliers(no_outliers,4)';

imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;

%% fiterror in patches in fov

Npatch = 20;
zrange = [-1000,1000];
%zrange = [0,1000];
%zrange = [-1000,0];


%% For the second fit, the range in z is smaller
%zrange = [-503.5 503.5];
zrange = [-503.5 0];
%zrange = [0 503.5];
%%
errM_fov = zeros(Npatch,Npatch);
errM_fov_beads = zeros(Npatch,Npatch);

for xi = 1:Npatch
    for yi = 1:Npatch

        xrange = [floor((imgSizeX/Npatch)*(xi-1))+1,floor((imgSizeX/Npatch))*xi];
        yrange = [floor((imgSizeY/Npatch)*(yi-1))+1,floor((imgSizeY/Npatch))*yi];
        indices = find((roixy(1,:) >= xrange(1) & roixy(1,:) <= xrange(2) & roixy(2,:) >= yrange(1) & roixy(2,:) <= yrange(2) & z_positions > zrange(1) & z_positions < zrange(2) ));
        %indices = indices(indices<1e5);
        Ncfg = numel(indices);

        % select spots with these indices
        allspots_selected = allspots(:,:,:,indices);
        
        allspots_modeled = mu_smlmaber;
        allspots_modeled_beads = mu_beadaber;

        allspots_modeled = allspots_modeled(:,:,:,indices);
        allspots_modeled_beads = allspots_modeled_beads(:,:,:,indices);

        errM = get_fiterror(allspots_modeled,allspots_selected,params);
        errM_beads = get_fiterror(allspots_modeled_beads,allspots_selected,params);

        % inspect the final result
        %subplot(Npatch,Npatch,(xi-1)*Npatch+yi);
        errM_fov(yi,xi) = mean(errM(1,:),"all");
        errM_fov_beads(yi,xi) = mean(errM_beads(1,:),"all");

        if(Ncfg == 0)
            errM_fov(yi,xi) = NaN;
            errM_fov_beads(yi,xi) = NaN;
        end
       
    end
end

%%
%imagesc(cat(2,errM_fov,errM_fov_beads))
figure 
xsize = 1e-3*params.imgSizeX*params.pixelsize;
ysize = 1e-3*params.imgSizeY*params.pixelsize;
xx = linspace(xsize/(2*Npatch),xsize - xsize/(2*Npatch),Npatch);
yy = linspace(ysize/(2*Npatch),ysize - ysize/(2*Npatch),Npatch);
fraction_chisquare = errM_fov_beads./errM_fov;
fraction_chisquare = errM_fov./errM_fov_beads;
img = imagesc(xx,yy,fraction_chisquare,[0.97,1.03]);
%img = imagesc(xx,yy,errM_fov_beads,[1,1.3]);%,[1,1.4]);
%img = imagesc(xx,yy,errM_fov,[1,1.3]);
%img = imagesc(xx,yy,errM_fov,[289,462]);%,[1,1.6]);
set(img,'AlphaData',~isnan(fraction_chisquare))
set(gca,'YDir','normal') 
xlabel('x (\mum)')
ylabel('y (\mum)')
axis equal
colormap(redblue)
%colormap parula
cbh = colorbar;

%%
% xlim([0 177])
% ylim([0 177])
%%
%set(cbh,'YTick',[0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05])
%set(cbh,'YTick',[0.90,0.95,1.00,1.05,1.10])
ax = gca;
ax.FontSize = 26;
% fontsize(16,'points')

mean(fraction_chisquare(~isnan(fraction_chisquare)),"all")
mean(errM_fov_beads(~isnan(errM_fov_beads)),"all")
mean(errM_fov(~isnan(errM_fov)),"all")



%% Calculate location arrows

data1 = data_initial_fit;
data2 = data_second_fit;

params = data1.params;
localizations_1 = data1.localizations;
localizations_2 = data2.localizations;
[C,index_1,index_2] = intersect(localizations_1(:,1),localizations_2(:,1));

xyz_1 = localizations_1(index_1,2:4);
xyz_2 = localizations_2(index_2,2:4);

%%
% Set parameters
imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;
pixelsize = params.pixelsize;

location_difference = xyz_2 - xyz_1;

%% define FOV patches
Npatch_xy = 15;
Npatch_z = 1;
minspots = 10; % Leave out patches with very few emitters

zmin = min(cat(1,localizations_1(:,4),localizations_2(:,4)),[],1);
zmax = max(cat(1,localizations_1(:,4),localizations_2(:,4)),[],1);

minxyz = [0 0 zrange(1)];
maxxyz = [imgSizeX*pixelsize imgSizeY*pixelsize zrange(2)];
Xpatch = minxyz(1)+(0:Npatch_xy)*(maxxyz(1)-minxyz(1))/Npatch_xy;
Ypatch = minxyz(2)+(0:Npatch_xy)*(maxxyz(2)-minxyz(2))/Npatch_xy;
Zpatch = minxyz(3)+(0:Npatch_z)*(maxxyz(3)-minxyz(3))/Npatch_z;

mean_location_difference = zeros(3,Npatch_xy,Npatch_xy,Npatch_z);
spots_per_patch = zeros(Npatch_xy,Npatch_xy,Npatch_z);

for ii = 1:Npatch_xy
    for jj = 1:Npatch_xy
        for kk=1:Npatch_z
            FOVfilter = (xyz_1(:,1)>Xpatch(ii))&(xyz_1(:,1)<Xpatch(ii+1))&...
                        (xyz_1(:,2)>Ypatch(jj))&(xyz_1(:,2)<Ypatch(jj+1))&...
                        (xyz_1(:,3)>Zpatch(kk))&(xyz_1(:,3)<Zpatch(kk+1));
    
            spots_per_patch(ii,jj,kk) = size(location_difference(FOVfilter,:),1);
            if spots_per_patch(ii,jj,kk) > minspots
                mean_location_difference(:,ii,jj,kk) = mean(location_difference(FOVfilter,:));
            % else
            %     mean_location_difference(:,ii,jj,kk) = NaN;
            end
        end
    end
end

%% Plot arrows (xy, 1 img for all z)

scalefac = 0.1; %6(2D); %0.2(3D);
x_legend = -5;
y_legend = -4;
u_legend = 100*scalefac; % k*scalefac means the legend arrow is k nm
v_legend = 0;

figure
x = 1e-3*Xpatch(1:end-1);
y = 1e-3*Ypatch(1:end-1);
[xx,yy] = meshgrid(x,y);
mean_location_difference_scaled = -1*scalefac*mean_location_difference;
%quiver3(xx(:)',yy(:)',zeros(size(xx(:)')),mean_location_difference_scaled(1,:),mean_location_difference_scaled(2,:),mean_location_difference_scaled(3,:),0,'LineWidth',2)
quiver(x_legend,y_legend,u_legend,v_legend,0,'LineWidth',3,'color',"#0072BD","MaxHeadSize",10)

hold on
quiver(xx(:)',yy(:)',mean_location_difference_scaled(2,:),mean_location_difference_scaled(1,:),0,'LineWidth',3,'AutoScale','off','color',"#0072BD","MaxHeadSize",10)
axis square
%%%%
x_pos = -x_legend-4;
y_pos = -y_legend-3;
text(x_pos , y_pos , ...
     '10 nm', 'Color', 'k', 'HorizontalAlignment', 'center', 'FontSize', 14);
fontsize(20,'points')
%%%%
%%
% xlim([-8 100]) % Lidke
% ylim([-8 100])
%%
%xlim([-8 177]) % Fu
%ylim([-8 177])
xlabel('x (\mum)')
ylabel('y (\mum)')
ax=gca;
ax.FontSize = 40;

% figure
% x = 1e-3*Xpatch(1:end-1);
% y = 1e-3*Ypatch(1:end-1);
% [xx,yy] = meshgrid(x,y);
% mean_location_difference_z_scaled = 0.4*mean_location_difference;
% quiver(xx(:)',yy(:)',mean_location_difference_z_scaled(1,:),mean_location_difference_z_scaled(2,:),0,'LineWidth',2)
% xlabel('x (\mum)')
% ylabel('y (\mum)')

%% Plot z difference as heat map over the FOV

figure 
xsize = 1e-3*params.imgSizeX*params.pixelsize;
ysize = 1e-3*params.imgSizeY*params.pixelsize;
xx = linspace(xsize/(2*Npatch),xsize - xsize/(2*Npatch),Npatch);
yy = linspace(ysize/(2*Npatch),ysize - ysize/(2*Npatch),Npatch);
img = imagesc(xx,yy,squeeze(mean_location_difference(3,:,:)),[-100 100]);
set(gca,'YDir','normal') 
xlabel('x (\mum)')
ylabel('y (\mum)')
axis equal
colormap(redblue)
%colormap parula
cbh = colorbar;
% xlim([0 177])
% ylim([0 177])

%%
% xlim([0 98])
% ylim([0 98])
%%
ax = gca;
ax.FontSize = 26;


for ii = 1: size(mean_location_difference,1)
    dim = mean_location_difference(ii,:,:);
    dim  = dim(~isnan(dim));
    disp(mean(dim,'all'))

end


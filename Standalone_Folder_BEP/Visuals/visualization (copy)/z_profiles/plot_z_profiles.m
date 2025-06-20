% Plot average z-profiles within a region

%% Load data
clear all
close all

path_input_data_1 = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\localizations_after_drift_correction\localizations_gain12_3okt.mat';
path_input_data_2 = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\localizations_after_drift_correction\localizations_gain12_constaber_3okt.mat';

data_1 = load(path_input_data_1);
data_2 = load(path_input_data_2);

x1 = data_1.localizations(:,2);
y1 = data_1.localizations(:,3);
z1 = data_1.localizations(:,4);

x2 = data_2.localizations(:,2);
y2 = data_2.localizations(:,3);
z2 = data_2.localizations(:,4);

params = data_1.params;
%% Select region
rangex = [0.0*params.pixelsize*params.imgSizeX 0.1*params.pixelsize*params.imgSizeX];
rangey = [0.72*params.pixelsize*params.imgSizeY 0.82*params.pixelsize*params.imgSizeY];

pixelsize = 5; % Should be ~half the size of the crlb. The smaller, the large the rendered image becomes.
alpha = 0.0;
beta = 0.0;
sigma = 2.0;
cmax = 1.0;
image_rendered = renderprojection(x1,y1,z1,alpha,beta,rangex,rangey,pixelsize,sigma,false,cmax);
imagesc(image_rendered,[0 cmax])
pbaspect([size(image_rendered,2) size(image_rendered,1) 1]);
colormap hot
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);

%%
figure
imshow(5*image_rendered)
h = impoly;
wait(h);
binaryMask = createMask(h);

%%
% First filter the coordinates that are in the region.
% For each coordinate, convert the x and y coordinates to a pixel coordinate in the mask
% If that coordinate is 1 in the binanymask, select the corresponding
% index.

region_indices = find(x1 >= rangex(1) & x1 <= rangex(2) & y1 >= rangey(1) & y1 <= rangey(2));
selected_indices = [1];
for idx=region_indices'
    xi = x1(idx);
    yi = y1(idx);

    xi_pixel = max(1,floor((xi-rangex(1))/pixelsize));
    yi_pixel = max(1,floor((yi-rangey(1))/pixelsize));

    if binaryMask(xi_pixel,yi_pixel)
        selected_indices = cat(1,selected_indices,idx);
    end
end
selected_indices = selected_indices(2:end);

%% Make z-plot
zrange = [-500 500];
zpixelsize = 10;
rz = zrange(1):zpixelsize:zrange(end);
zbin_center = zrange(1)+zpixelsize/2:zpixelsize:zrange(end)-zpixelsize/2;

z_selected_1 = z1(selected_indices);
z_selected_2 = z2(selected_indices);

z_hist_1 = histcounts(z_selected_1,rz);
z_hist_2 = histcounts(z_selected_2,rz);

figure
plot(zbin_center,z_hist_1,'LineWidth',3)
hold on
plot(zbin_center,z_hist_2,'LineWidth',3)
xlabel('z(nm)')
ylabel('Counts')
legend('Fitted aberrations','Constant astigmatism')
grid on
ax=gca;
ax.FontSize = 20;



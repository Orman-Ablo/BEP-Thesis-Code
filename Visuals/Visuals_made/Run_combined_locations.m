% input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Combining Results\combined_locations_second_fit.mat";
% input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Combining Results\combined_locations_initial_fit.mat";
% load(input,'locations')

input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\data\Combined\localization_combined_second_fit_scaled.mat";
%input = "C:\Users\Naam\Documents\BEP\data\Combined\localization_combined_initial_fit.mat";
params = set_parameters_microtubules_3D;

load(input,'localizations_drift_corrected')
locations = localizations_drift_corrected(:,2:4);

ncolors = 19;
colorrange = [-500,400];

xlim = [0.0 1.0];
ylim = [0.0 1.0];
pixelsize = 50
locations(:,3) = -locations(:,3);
sigma = 2.0;
clim = 0.95; % the brightest (1-clim) fraction of pixels gets clipped 
figure;
img = render_3D_func(locations,ncolors,colorrange,xlim,ylim,pixelsize,sigma,clim,params);

xlabel('z-position (nm)')
%%%%%
% img = rgb2gray(img);
% imshow(img);
%%%%%

% my_colormap = jet;
% colormap(my_colormap);
% cb = colorbar('southoutside');
% cb.Ticks = [-500,-400,-300,-200,-100,0,100,200,300,400];
scale_bar = 1e4; % 10 micrometers


scale_bar_length_px = round((scale_bar) / pixelsize);
hold on
y_pos = size(img, 1) - 20; % 20 pixels from bottom
x_pos = 0 + scale_bar_length_px ; % 20 pixels from right

% Draw a white rectangle as a scalebar (1 pixel height)
rectangle('Position', [x_pos, y_pos, scale_bar_length_px, 10], ...
          'FaceColor', 'w', 'EdgeColor', 'none');

% Optional: Add text label
text(x_pos + scale_bar_length_px/2, y_pos -40, ...
     '10 μm', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', 22);
ax = gca;
ax.FontSize = 26;





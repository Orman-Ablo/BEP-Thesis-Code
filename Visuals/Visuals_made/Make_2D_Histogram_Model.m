input = 'C:\Users\Naam\Documents\BEP\BEP Code Roman\input data\transfer_3055388_files_4e4cf1e2\localization_segmentation_roi17_thr15p5_Data0001.mat';
load(input,'allspots','localizations_with_outliers','mu')
% input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\data\final\localization_OTF_localization_segmentation_roi17_thr15p5_Data0001.mat";
% load(input,'Spots_fit','localizations_with_outliers','mu')
% allspots = Spots_fit;

jcfg = 10;

data = (allspots(:,:,1,jcfg)-localizations_with_outliers(jcfg,end))/localizations_with_outliers(jcfg,end-1);

% %% Produce a gaussian model
% sigma = 2.3;
% A = 2.5;
% x= linspace (-10,10,17);
% y= linspace(-10,10,17);
% [X, Y] = meshgrid(x, y);
% gauss = A/(2*pi*sigma^2)*exp(-(X.^2+Y.^2)/(2*sigma^2));
% size(gauss)






model = (mu(:,:,1,jcfg)-localizations_with_outliers(jcfg,end))/localizations_with_outliers(jcfg,end-1);

% Create a 3D bar plot
figure;


% either use params.pixelsize, or change the pixelsize outright
pixelsize = 95;

% Mx = size(data,1);
% My = size(data,2);
% 
% xlims = [-pixelsize*Mx/2 pixelsize*Mx/2];
% ylims = [-pixelsize*My/2 pixelsize*My/2];
% 
% x = linspace(-pixelsize*Mx/2, pixelsize*Mx/2, 17);
% y = linspace(-pixelsize*My/2, pixelsize*My/2, 17);

% colormap over intensity instead of position
h1 = bar3(data,1);
for k = 1:length(h1)
    % xdata = zeros(size(zdata));
    % ydata = zeros(size(zdata));
    % 
    % % Create meshgrid for x and y coordinates
    % [xdata, ydata] = meshgrid(x(k), y);
    % 
    % % Update XData and YData
    % h(k).XData = xdata;
    % h(k).YData = ydata;

    % Set the color data for each bar
    zdata = h1(k).ZData;
    h1(k).CData = zdata;
    h1(k).FaceColor = 'interp';
end
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Bar plot of Data Spot');

figure;
h2 =bar3(model,1);
for k = 1:length(h2)
    zdata = h2(k).ZData;
    h2(k).CData = zdata;
    h2(k).FaceColor = 'interp';
end
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Bar Plot of Model Spot');
function [rendered_image] = render_3D_func(locations,ncolors,colorrange,xlim,ylim,pixelsize,sigma,clim,params)
% Creates a rendered image that is color coded for z
% The color is given by the weighted average depth at that x,y pixel
% location
% Author: Isabel Droste, TU Delft, 2024

    Ncfg = size(locations,1);
    
    theta_local_zset = zeros(Ncfg,4);
    theta_local_zset(:,1:3) = locations;
    
    dcolor = (colorrange(2) - colorrange(1))/(ncolors-1);
    theta_local_zset(:,4) = min(max(floor((theta_local_zset(:,3)-colorrange(1))/dcolor),0),ncolors-1);
    
    dx = 0;
    dy = 0;
    
    % corrections for fd deeploc
    %dx = -0.0168; %fig 3cde
    %dy = -0.0168;

    %dx = -0.011; %sub fig 1bcd
    %dy = -0.0107;

    %dx = -0.009; %sub fig 2bcd
    %dy = -0.0025;

    rangex = [(xlim(1)+dx)*params.pixelsize*params.imgSizeX (xlim(2)+dx)*params.pixelsize*params.imgSizeX];
    rangey = [(ylim(1)+dy)*params.pixelsize*params.imgSizeY (ylim(2)+dy)*params.pixelsize*params.imgSizeY];
    
    npix = size(rangex(1):pixelsize:rangex(end),2)-1; % calculate automatically
    npiy = size(rangey(1):pixelsize:rangey(end),2)-1;
    final_image = zeros(npix,npiy,ncolors);

    % Create image for each color level. Then calculate pixelwise average
    % color.
    for n=1:ncolors
        %n
        selected_indices = (theta_local_zset(:,4) == n-1);
    
        locations_selected = locations(selected_indices,:);
        %size(locations_selected)
        fprintf('Processing depth %i\n',n)
        y = locations_selected(:,1);
        x = locations_selected(:,2);
        z = locations_selected(:,3);
    
        alpha = 0.0;
        beta = 0.0;
        cmax = 1;
        image_rendered = renderprojection(x,y,z,alpha,beta,rangex,rangey,pixelsize,sigma,false,cmax);
        final_image(:,:,n) = image_rendered;
    end
    
    count_up = zeros(1,1,ncolors);
    count_up(1:ncolors) = 0:ncolors-1;
    count_up_matrix = repmat(count_up,npix,npiy);
    
    multiplied_final = (final_image.*count_up_matrix)./sum(final_image,3);
    color_weights = sum(multiplied_final,3);
    color_weights(isnan(color_weights)) = 0;
    
    nc = floor(256*6/5);
    cmap = hsv(nc);
    my_colormap = colormap(cmap(1:ceil(nc*5/6),:));
    mgint = floor(color_weights*255/(ncolors-1));
    rgb = ind2rgb(mgint,my_colormap);

    % clip top part using clim
    rescaled_img = rescale(sum(final_image,3));
    pct = prctile(rescaled_img(:),clim*100);
    clipped_img = min(rescaled_img/pct,1);
    rendered_image = permute(rgb.*clipped_img,[2 1 3]);
    

    
    
    imagesc(rendered_image)
    pbaspect([size(rendered_image,2) size(rendered_image,1) 1]);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);

    colormap(my_colormap);
    caxis([-500 400]); 
    cb = colorbar('southoutside');
    cb.TickLabels = [-500,-200,0,200,400];
end


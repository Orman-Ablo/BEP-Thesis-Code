function [rendered_image] = render_3D_crossection_func(locations,colorrange,xlim,ylim,pixelsize,sigma,clim,params)
% Creates a rendered cross section that is color coded for z
% Author: Isabel Droste, TU Delft, 2024

    dx = 0;
    dy = 0;
    
    % corrections for fd deeploc
    %dx = -0.0168; %3c
    %dy = -0.0168;

    %dx = -0.011; %sub fig 1bcd
    %dy = -0.0107;

    %dx = -0.009; %sub fig 2bcd
    %dy = -0.0025;
    
    % select spots in x-range
    selected_indices = (locations(:,1) >= (ylim(1)+dy)*params.pixelsize*params.imgSizeY & locations(:,1) <= (ylim(2)+dy)*params.pixelsize*params.imgSizeX);
    sum(selected_indices)
    xx = locations(selected_indices,1);
    yy = locations(selected_indices,2);
    zz = locations(selected_indices,3);
    
    % Create 2D image
    rangeyy = [(xlim(1)+dx)*params.pixelsize*params.imgSizeY (xlim(2)+dx)*params.pixelsize*params.imgSizeY];
    rangezz = colorrange;
    alpha = 0.0;
    beta = 0.0;
    cmax = 1.0;
    img = renderprojection(yy,zz,xx,alpha,beta,rangeyy,rangezz,pixelsize,sigma,false,cmax);
    
    % add color for z in xz plot
    nc = floor(256*6/5);
    cmap = hsv(nc);
    my_colormap = colormap(cmap(1:ceil(nc*5/6),:));
    
    size_y = size(img,1);
    size_z = size(img,2);
    colorgrad = zeros(size_y,size_z,3);
    for i=1:size_z
        idx =  max(floor(i*size(my_colormap,1)/size_z),1);
        colorgrad(:,i,:) = repmat(my_colormap(idx,:),[size_y,1]);
    end
    
    rescaled_img = rescale(img);
    pct = prctile(rescaled_img(:),clim*100)
    clipped_img = min(rescaled_img/pct,1);
    rendered_image = colorgrad.*clipped_img;

    %rendered_image = colorgrad.*repmat(img,[1,1,3]);
    imagesc(rangeyy,rangezz,permute(rendered_image,[2 1 3]));
    pbaspect([size(rendered_image,1) size(rendered_image,2) 1]);
    xlabel('nm')
    ylabel('nm')
    set(gca,'YDir','normal')
end
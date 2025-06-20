function [fov_coordinates_normalized,fov_coordinates_physical] = get_fov_coordinates(roixy,x,y,params)
%get_fov_coordinates Returns normalized and physical FOV coordinates. The
%normalized coordinates are in the square interval [-1,1] x [-1,1] and the
%physical coordinates are in um.
pixelsize = params.pixelsize;
xsize = params.pixelsize*params.imgSizeX; % in pysical coordinates
ysize = params.pixelsize*params.imgSizeY;

fov_coordinates_physical_x = (roixy(1,:)*pixelsize + x)*1e-03;
fov_coordinates_physical_y = (roixy(2,:)*pixelsize + y)*1e-03;
fov_coordinates_physical = cat(1,fov_coordinates_physical_x,fov_coordinates_physical_y);

fov_coordinates_x = -1+2*(roixy(1,:)*pixelsize + x)/xsize;
fov_coordinates_y = -1+2*(roixy(2,:)*pixelsize + y)/ysize;
fov_coordinates_normalized = cat(1,fov_coordinates_x,fov_coordinates_y);

end
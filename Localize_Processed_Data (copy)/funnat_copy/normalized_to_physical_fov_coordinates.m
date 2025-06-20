function [Xim,Yim] = normalized_to_physical_fov_coordinates(Xn,Yn,params);
%normalized_to_physical_fov_coordinates Converts normalized fov coordinates
%in the interval [-1,1] to the physical FOV size.

pixelsize = params.pixelsize;
imgSizeX = params.imgSizeX;
imgSizeY = params.imgSizeY;

Xim = 1e-3*(Xn+1)*pixelsize*imgSizeX/2;
Yim = 1e-3*(Yn+1)*pixelsize*imgSizeY/2;

end
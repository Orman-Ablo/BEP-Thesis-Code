function PSFsum = sum_shift_ups_PSFs(allspots,allpositions,upsfac)
% This function is for upsampling, shifting to the fitted spot centers, and
% subsequently summing a stack of measured single molecule spots. Only
% integer upsampling is allowed and we use zero padding in Fourier space to
% achieve this.
%
% copyright Sjoerd Stallinga, TU Delft, 2023

[Nx,Ny,Nstack] = size(allspots);
Mx = upsfac*Nx;
My = upsfac*Ny;
centerNx = floor(Nx/2)+1;
centerNy = floor(Ny/2)+1;
centerMx = floor(Mx/2)+1;
centerMy = floor(My/2)+1;
xrange = centerMx-centerNx+1:centerMx-centerNx+Nx;
yrange = centerMy-centerNy+1:centerMy-centerNy+Ny;

% spatial frequency vectors, this only works if the set of positions is
% measured in pixel units
qxx = ((1:Nx)-centerNx)/Nx;
qyy = ((1:Ny)-centerNy)/Ny;
[Qx,Qy] = meshgrid(qyy,qxx);

% loop over stack to upsample, shift, and sum
PSFsum = zeros(upsfac*Nx,upsfac*Ny);
for js = 1:Nstack
  spot = allspots(:,:,js);
  ft_spot = fftshift(fft2(spot));
  shiftpos = allpositions(:,js);
  shiftpos = shiftpos-(upsfac-1)/2/upsfac;
  phase = 2*pi*(Qx*shiftpos(1)+Qy*shiftpos(2));
  ft_spot = exp(1i*phase).*ft_spot;
  ft_spot_extend = zeros(Mx,My);
  ft_spot_extend(xrange,yrange) = ft_spot;
  ups_spot = real(ifft2(ifftshift(ft_spot_extend)));
%   ups_spot = real(ifft2(ft_spot_extend));
  PSFsum = PSFsum+ups_spot;
end

end
function PSFsum = sum_shift_ups_PSFs_mod(allspots,allpositions,upsfac,bg,Nph)
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
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Method 1: General scaling of the background
  bg_spot = bg(js)/upsfac^2; % in order to correspond to the change through interpolation, the background needs to be scaled DOWN
  Nph_spot = Nph(js);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % %Method 2: recomputing the background
  % 
  % Using the method in initialvalues_phasor? (more robust)
  % % initial estimate of background bg from median value of the rim pixels,
  % % and signal photon count Nph from total photon count in ROI
  %   rimpixels = zeros(2*Mx+2*My-4,1);
  %   rimpixels(1:Mx-1) = fxy(1:Mx-1,1); % fxy = allspots
  %   rimpixels(Mx:2*Mx-2) = fxy(2:Mx,My);
  %   rimpixels(2*Mx-1:2*Mx+My-3) = fxy(Mx,1:My-1);
  %   rimpixels(2*Mx+My-2:2*Mx+2*My-4) = fxy(1,2:My);
  %   bg0 = median(rimpixels);
  %   bg0 = max(bg0,1);
  %   Nph0 = sum(sum(fxy))-Mx*My*bg0;
  %   if (Nph0<0)
  %       Nph0 = sum(sum(fxy));
  %   end

  % % refinement of initial estimate of photon count and background by linear
  %   % regression, i.e. a linear least squares fit to n = Nph*PSF+bg
  %   % this procedure is entirely free from ad-hoc parameters and gives a bias
  %   % free initial estimate
  % 
  %   % first we compute the PSF for this initial position estimate (x0,y0)
  %   params_tmp = params;
  %   params_tmp.xemit = x0;
  %   params_tmp.yemit = y0;
  %   params_tmp.zemit = z0;
  %   params_tmp.K = 1;
  %   [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = ...
  %     get_pupil_matrix(params_tmp,zernikeCoefficients(jcfg,:)); %%%% new : extra input value
  %   [FieldMatrix,FieldMatrixDerivatives] = ...
  %     get_field_matrix_derivatives(params_tmp,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
  %   [PSF,~] = get_psfs_derivatives(params_tmp,PupilMatrix,FieldMatrix,FieldMatrixDerivatives);
  % 
  %   % use this for a refinement with weighted linear least squares
  %   % which is here mathematically equivalent to solving the MLE equations
  %   % d(log(L))/dNph=0 and d(log(L))/dbg=0 for Nph and bg
  %   % this gives an additional improvement in number of MLE iterations
  %   weight = 1./(Nph0*PSF+bg0);
  %   % weight = ones(size(PSF));
  %   Hav = mean(weight(:).*PSF(:));
  %   H2av = mean(weight(:).*PSF(:).^2);
  %   nav = mean(weight(:).*fxy(:));
  %   nHav = mean(weight(:).*fxy(:).*PSF(:));
  %   wav = mean(weight(:));
  %   detty = H2av*wav-Hav^2;
  %   Nph = (nHav*wav-nav*Hav)/detty;
  %   bg = (nav*H2av-nHav*Hav)/detty;
  % 
  %   % provision for negative or too large estimates
  %   if (bg < 1 || bg > 5*bg0)
  %       Nph = Nph+Mx*My*(bg-bg0);
  %       bg = bg0;
  %   end
  %   if (Nph < 0 || Nph > 2*sum(sum(fxy)))
  %       Nph = Nph0;
  %   end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ft_spot = fftshift(fft2(spot));
  shiftpos = allpositions(:,js);
  shiftpos = shiftpos-(upsfac-1)/2/upsfac;
  phase = 2*pi*(Qx*shiftpos(1)+Qy*shiftpos(2));
  ft_spot = exp(1i*phase).*ft_spot;
  ft_spot_extend = zeros(Mx,My);
  ft_spot_extend(xrange,yrange) = ft_spot;
  ups_spot = real(ifft2(ifftshift(ft_spot_extend)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ups_spot = (ups_spot - bg_spot)/Nph_spot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%   ups_spot = real(ifft2(ft_spot_extend));
  PSFsum = PSFsum+ups_spot;
end

end
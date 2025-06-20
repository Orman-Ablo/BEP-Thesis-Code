function OTF = get_otf(PSF,parameters)
% This function calculates the 2D/3D OTF from the 2D/3D PSF by CZT.
%
% copyright Sjoerd Stallinga, TU Delft, 2021-2023

NA = parameters.NA;
lambda = parameters.lambda;
refin = parameters.refimmnom; % take nominal coverslip/immersion fluid refractive index
Notf = parameters.Notf;
Notfz = parameters.Notfz;
Mpsfx = parameters.Mpsfx;
Mpsfy = parameters.Mpsfy;
xrange = parameters.psfxrange;
yrange = parameters.psfyrange;

% PSF and OTF space size (in diffraction units)
PSFSizex = xrange*NA/lambda;
PSFSizey = yrange*NA/lambda;
%%
% Rescaling PSFsize to match the size of the ROIs
PSFSizex = (parameters.Mx*parameters.pixelsize/2)*NA/lambda;% This should match the Sizes of the ROIs
PSFSizey = (parameters.My*parameters.pixelsize/2)*NA/lambda;%
%%
OTFSize = 2.0; % OTF cutoff in xy-direction
switch parameters.fitmodel
  case 'xyz'
    
    Mpsfz = parameters.Mpsfz;
    % incorrect, arbitrary z-range
    PSFSizez = parameters.psfzrange*(refin-sqrt(refin^2-NA^2))/lambda*2;
    % correct z range
    PSFSizez = parameters.z_range*(refin-sqrt(refin^2-NA^2))/lambda; % this is the upper bound of the Z-window.
    % larger range than default cutoff (size=1) to accommodate possible RI mismatch
    OTFSizez = 2.0; 
end

% % apply Tukey windowing to computed PSF to prevent ringing effects
% cosfrac = 0.25;
% winx = tukeywin(Mpsfx,cosfrac);
% winy = tukeywin(Mpsfy,cosfrac);
% switch parameters.fitmodel
%   case 'xy'
%     [winY,winX] = meshgrid(winy,winx);
%     totwinXY = winX.*winY;
%     PSF = totwinXY.*PSF;
%   case 'xyz'
%     winz = tukeywin(Mpsfz,cosfrac);
%     [winY,winX,winZ] = meshgrid(winy,winx,winz);
%     totwinXYZ = winX.*winY.*winZ;
%     PSF = totwinXYZ.*PSF;
% end

% calculate auxiliary vectors for chirpz
switch parameters.fitmodel
  case 'xy'
%     [Ax,Bx,Dx] = prechirpz(PSFSizex,OTFSize,Mpsfx,Notf);
%     [Ay,By,Dy] = prechirpz(PSFSizey,OTFSize,Mpsfy,Notf);
    allPSFsize = [PSFSizex PSFSizey];
    allOTFSize = [OTFSize OTFSize];
    allMpsf = [Mpsfx Mpsfy];
    allNotf = [Notf Notf];
  case 'xyz'
    allPSFsize = [PSFSizex PSFSizey PSFSizez];
    allOTFSize = [OTFSize OTFSize OTFSizez];
    allMpsf = [Mpsfx Mpsfy Mpsfz];
    allNotf = [Notf Notf Notfz];
end
[allA,allB,allD] = prechirpzn(allPSFsize,allOTFSize,allMpsf,allNotf);

% calculation of OTF by CZT from PSF
% IntermediateImage = transpose(cztfunc(PSF,Ay,By,Dy));
% OTF = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));

OTF = cztfuncn(PSF,allA,allB,allD,true);

% OTF = cztfuncn2(PSF,allA,allB,allD,allPSFsize,allOTFSize,allMpsf,allNotf,true);

OTF = conj(OTF); % needed to make FT and IFT consistent
OTF = flip(OTF,3);

% additional normalization of OTF not needed!

% % apply Tukey windowing to computed OTF to prevent ringing effects
% cosfrac = 0.25;
% winx = tukeywin(Notf,cosfrac);
% winy = tukeywin(Notf,cosfrac);
% switch parameters.fitmodel
%   case 'xy'
%     [winY,winX] = meshgrid(winy,winx);
%     totwinXY = winX.*winY;
%     OTF = totwinXY.*OTF;
%   case 'xyz'
%     winz = tukeywin(Notfz,cosfrac);
%     [winY,winX,winZ] = meshgrid(winy,winx,winz);
%     totwinXYZ = winX.*winY.*winZ;
%     OTF = totwinXYZ.*OTF;
% end

end



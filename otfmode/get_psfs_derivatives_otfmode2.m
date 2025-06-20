function [PSF,PSFderivatives] = get_psfs_derivatives_otfmode2(OTF,parameters,compders)
% This function calculates the PSF and the derivatives w.r.t. the (x,y) 
% and possibly z, position of the emitter using FT from a pre-computed OTF,
% which must be passed as an argument (Notf x Notf array for fitmodel = 'xy'
% or Notf x Notf x Notfz for fitmodel = 'xyz') to the function.
% So far, it only applies to fitting a single 2D ROI, e.g. multi-focal
% fitting on a set of 2D ROIs is not implemented (Mz = 1).
%
% copyright Sjoerd Stallinga, TU Delft, 2021-2023
%
% flag compders indicates mode with/without computation PSFderivatives 



%% same as 'get_psfs_derivatives_otfmode', copy function

xemit = parameters.xemit;
yemit = parameters.yemit;
NA = parameters.NA;
lambda = parameters.lambda;
refin = parameters.refimmnom;
Notf = parameters.Notf;
Mx = parameters.Mroix;
My = parameters.Mroiy;
xrange = parameters.xroirange;
yrange = parameters.yroirange;
switch parameters.fitmodel
  case 'xy'
    numders = 2;
  case 'xyz'
    Notfz = parameters.Notfz;
    numders = 3;
    switch parameters.ztype
      case 'medium'
        zemit = parameters.zemit;
      case 'stage'
        zemit = parameters.zstage;
    end
end

% OTF space size and sampling (in diffraction units)
OTFSize = 2.0; % OTF cutoff in xy-direction
DxyOTF = 2*OTFSize/Notf;
XYOTF = -OTFSize+DxyOTF/2:DxyOTF:OTFSize;
[YOTF,XOTF] = meshgrid(XYOTF,XYOTF);

% calculate auxiliary vectors for chirpz
% [Ax,Bx,Dx] = prechirpz(OTFSize,ROISizex,Notf,Mx);
% [Ay,By,Dy] = prechirpz(OTFSize,ROISizey,Notf,My);
allA = parameters.allA;
allB = parameters.allB;
allD = parameters.allD;

% compute 2D-OTF from 3D-OTF and z-derivative in case fitmodel='xyz'
switch parameters.fitmodel
  case 'xyz'
    OTFSizez = 2.0; % larger range than default cutoff (size=1) to accommodate possible RI mismatch
    DzOTF = 2*OTFSizez/Notfz;
    
    ZOTF = -OTFSizez+DzOTF/2:DzOTF:OTFSizez;

    wavevectorz = (2*pi*(refin-sqrt(refin^2-NA^2))/lambda)*ZOTF;
    ZPositionPhaseMask = exp(-1i*zemit*wavevectorz);
    ZPositionPhaseMaskrep = repmat(ZPositionPhaseMask,Notf^2,1);
    ZPositionPhaseMaskrep = reshape(ZPositionPhaseMaskrep,[Notf Notf Notfz]);
    OTFtemp = DzOTF*sum(OTF.*ZPositionPhaseMaskrep,3);


    if compders
      diffZPositionPhaseMask = -1i*wavevectorz.*exp(-1i*zemit*wavevectorz);
      diffZPositionPhaseMaskrep = repmat(diffZPositionPhaseMask,Notf^2,1);
      diffZPositionPhaseMaskrep = reshape(diffZPositionPhaseMaskrep,[Notf Notf Notfz]);
      zOTF = DzOTF*sum(OTF.*diffZPositionPhaseMaskrep,3);

    end
    OTF = OTFtemp;

    

end

% calculation of wavevector components in x and y and of emitter position
% dependent phase factor
wavevectorxy = zeros(Notf,Notf,2);
wavevectorxy(:,:,1) = (2*pi*NA/lambda)*XOTF;
wavevectorxy(:,:,2) = (2*pi*NA/lambda)*YOTF;
Wpos = xemit*wavevectorxy(:,:,1)+yemit*wavevectorxy(:,:,2);
PositionPhaseMask = exp(-1i*Wpos);
OTF = PositionPhaseMask.*OTF;


% calculation of PSF by CZT from OTF
% IntermediateImage = transpose(cztfunc(OTF,Ay,By,Dy));
% PSF = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
PSF = cztfuncn(OTF,allA,allB,allD,false);

%
% PSF = ifft2( fftshift( OTF ) );  
% PSF = ifftshift(PSF);
%
PSF = real(PSF);




% calculation of PSF derivatives, if so required
if compders
  PSFderivatives = zeros(Mx,My,1,numders);
  for jder = 1:2
    PupilFunction = -1i*wavevectorxy(:,:,jder).*OTF;
%     IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
%     PSFderivatives(:,:,1,jder) = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
    PSFderivatives(:,:,1,jder) = cztfuncn(PupilFunction,allA,allB,allD,false);
  end
  switch parameters.fitmodel
    case 'xyz'
      if compders
        PupilFunction = PositionPhaseMask.*zOTF;
%         IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
%         PSFderivatives(:,:,1,3) = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
        PSFderivatives(:,:,1,3) = cztfuncn(PupilFunction,allA,allB,allD,false);



      end
  end
  PSFderivatives = real(PSFderivatives);
  

else
  PSFderivatives = 0; % dummy value
end

end


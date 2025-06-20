function [XImage,YImage,FieldMatrix,FieldMatrixDerivatives] = ...
  get_field_matrix_derivatives(PupilMatrix,wavevector,wavevectorzimm,Waberration,allzernikes,parameters,compders)
% This function calculates the field matrix A_{jk}, which gives the j-th
% electric field component proportional to the k-th dipole vector
% component, as well as the derivatives of A_{jk} w.r.t. the xyz coordinates
% of the emitter and w.r.t. the emission wavelength lambda.
%
% copyright Sjoerd Stallinga, TU Delft, 2017
%
% parameters: NA, refractive indices of medium, wavelength (in nm), 
% nominal emitter position (in nm) with z-position from
% cover slip-medium interface, spot footprint (in nm), axial range (in nm),
% sampling in pupil with (even), sampling in image plane (odd), sampling in
% axial direction
%
% flag compders indicates mode with/without computation FieldMatrixDerivatives 

NA = parameters.NA;
lambda = parameters.lambda;
refmed = parameters.refmed;
xemit = parameters.xemit;
yemit = parameters.yemit;
zemit = parameters.zemit;
Npupil = parameters.Npupil;
switch parameters.fitmode
  case 'OTF'
    xrange = parameters.psfxrange;
    yrange = parameters.psfyrange;
    zmin = -parameters.psfzrange;
    zmax = parameters.psfzrange;
    zstage = 0;
    Mx = parameters.Mpsfx;
    My = parameters.Mpsfy;
    Mz = parameters.Mpsfz;
  case 'PSF'
    xrange = parameters.xrange;
    yrange = parameters.yrange;
    zmin = parameters.zrange(1);
    zmax = parameters.zrange(2);
    zstage = parameters.zstage;
    Mx = parameters.Mx;
    My = parameters.My;
    Mz = parameters.Mz;
end

% pupil and image size (in diffraction units)
PupilSize = 1.0;
ImageSizex = xrange*NA/lambda;
ImageSizey = yrange*NA/lambda;
ImageSizez = (zmax-zmin)/2;

% image coordinate sampling (in physical length units).
DxImage = 2*ImageSizex/Mx;
DyImage = 2*ImageSizey/My;
ximagelin = -ImageSizex+DxImage/2:DxImage:ImageSizex;
yimagelin = -ImageSizey+DyImage/2:DyImage:ImageSizey;
[YImage,XImage] = meshgrid(yimagelin,ximagelin);
XImage = (lambda/NA)*XImage;
YImage = (lambda/NA)*YImage;

% enlarging ROI in order to accommodate convolution with bead in
% computation PSF and PSFderivatives
if isfield(parameters,'bead')
  if parameters.bead == true
    beaddiameter = parameters.beaddiameter*NA/lambda;
    DeltaMx = 2*ceil(beaddiameter/DxImage);
    Mx = Mx+DeltaMx;
    ImageSizex = ImageSizex+DeltaMx*DxImage/2;
    DeltaMy = 2*ceil(beaddiameter/DyImage);
    My = My+DeltaMy;
    ImageSizey = ImageSizey+DeltaMy*DyImage/2;
    ximagelin = -ImageSizex+DxImage/2:DxImage:ImageSizex;
    yimagelin = -ImageSizey+DyImage/2:DyImage:ImageSizey;
    [YImage_tmp,XImage_tmp] = meshgrid(yimagelin,ximagelin);
    XImage_tmp = (lambda/NA)*XImage_tmp;
    YImage_tmp = (lambda/NA)*YImage_tmp;
  else
    XImage_tmp = XImage;
    YImage_tmp = YImage;
  end
else
  XImage_tmp = XImage;
  YImage_tmp = YImage;
end

% calculate auxiliary vectors for chirpz
% [Ax,Bx,Dx] = prechirpz(PupilSize,ImageSizex,Npupil,Mx);
% [Ay,By,Dy] = prechirpz(PupilSize,ImageSizey,Npupil,My);
allPupilSize = [PupilSize PupilSize];
allImageSize = [ImageSizex ImageSizey];
allNpupil = [Npupil Npupil];
allM = [Mx My];
[allA,allB,allD] = prechirpzn(allPupilSize,allImageSize,allNpupil,allM);

% calculation Zernike mode normalization
if strcmp(parameters.fitmodel,'aberrations')
  orders = parameters.aberrations(:,1:2);
  normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));
end

% set number of relevant parameter derivatives
switch parameters.fitmodel
  case 'xy'
    numders = 2;
  case 'xyz'
    numders = 3;
  case 'xylambda'
    numders = 3;
  case 'xyzlambda'
    numders = 4;
  case 'aberrations'
    numzers = size(parameters.aberrations,1);
    numders = 3+numzers;
end

% calculate nominal intensity normalization by flow of energy through lens
% aperture in case of SAF
if NA>refmed
  [normint_free,normint_fixed] = get_normalization(PupilMatrix,parameters);
  switch parameters.dipoletype
    case 'free'
      normfac0 = normint_free;
    case 'fixed'
      normfac0 = normint_fixed;
  end
end

% loop over emitter or stage z-position
if Mz==1
  ZImage = (zmin+zmax)/2;
else
  DzImage = 2*ImageSizez/Mz;
% enlarging axial range in order to accommodate convolution with bead
  if isfield(parameters,'bead')
    if parameters.bead == true
    beaddiameter = parameters.beaddiameter;
    DeltaMz = 2*ceil(beaddiameter/DzImage);
    Mz = Mz+DeltaMz;
    zmin = zmin-DeltaMz*DzImage/2;
    zmax = zmax+DeltaMz*DzImage/2;
    end
  end
  ZImage = zmin+DzImage/2:DzImage:zmax;
end
FieldMatrix = cell(2,3,Mz);
FieldMatrixDerivatives = cell(2,3,Mz,numders);

for jz = 1:numel(ZImage)
  zemitrun = ZImage(jz);  

% phase contribution due to position of the emitter
  Wpos = xemit*wavevector{1}+yemit*wavevector{2}+zemit*wavevector{3};
  Wpos = Wpos+zstage*wavevectorzimm;
  if strcmp(parameters.ztype,'medium')
    Wpos = Wpos+zemitrun*wavevector{3};
  end
  if strcmp(parameters.ztype,'stage')
    Wpos = Wpos+zemitrun*wavevectorzimm;
  end
  PositionPhaseMask = exp(-1i*Wpos);
  
 % needs extra normalization here per z-slice in case of SAF NA>refmed
  if NA>refmed
    PupilMatrix_temp = cell(2,3);
    for itel = 1:2
      for jtel = 1:3
        PupilMatrix_temp{itel,jtel} = PositionPhaseMask.*PupilMatrix{itel,jtel};
      end
    end
    [normint_free,normint_fixed] = get_normalization(PupilMatrix_temp,parameters);
    switch parameters.dipoletype
      case 'free'
        normfac = normint_free;
      case 'fixed'
        normfac = normint_fixed;
    end
    for itel = 1:2
      for jtel = 1:3
        PupilMatrix{itel,jtel} = sqrt(normfac0/normfac).*PupilMatrix{itel,jtel};
      end
    end
  end
  
  for itel = 1:2
    for jtel = 1:3
      % pupil functions and FT to matrix elements
      PupilFunction = PositionPhaseMask.*PupilMatrix{itel,jtel};
%       IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
%       FieldMatrix{itel,jtel,jz} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
      FieldMatrix{itel,jtel,jz} = cztfuncn(PupilFunction,allA,allB,allD,false);
      
      if compders
        % pupil functions for xy-derivatives and FT to matrix elements
        for jder = 1:2
          PupilFunction = -1i*wavevector{jder}.*PositionPhaseMask.*PupilMatrix{itel,jtel};
%           IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
%           FieldMatrixDerivatives{itel,jtel,jz,jder} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
          FieldMatrixDerivatives{itel,jtel,jz,jder} = cztfuncn(PupilFunction,allA,allB,allD,false);
        end

        % pupil functions for z-derivative and FT to matrix elements
        switch parameters.fitmodel
          case {'xyz','xyzlambda'}
            zderindex = 3;
            switch parameters.ztype
              case 'stage'
                PupilFunction = -1i*wavevectorzimm.*PositionPhaseMask.*PupilMatrix{itel,jtel};
              case 'medium'
                PupilFunction = -1i*wavevector{zderindex}.*PositionPhaseMask.*PupilMatrix{itel,jtel};
            end
%             IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
%             FieldMatrixDerivatives{itel,jtel,jz,zderindex} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
            FieldMatrixDerivatives{itel,jtel,jz,zderindex} = cztfuncn(PupilFunction,allA,allB,allD,false);
        end

        % pupil functions for lambda-derivative and FT to matrix elements
        switch parameters.fitmodel
          case 'xylambda'
            lambdaderindex = 3;
          case 'xyzlambda'
            lambdaderindex = 4;
        end
        
        switch parameters.fitmodel
          case {'xylambda','xyzlambda'}
            PupilFunction = -(2*pi*1i*Waberration/lambda^2).*PositionPhaseMask.*PupilMatrix{itel,jtel};
%             IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
%             FieldMatrixDerivatives{itel,jtel,jz,lambdaderindex} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
            FieldMatrixDerivatives{itel,jtel,jz,lambdaderindex} = cztfuncn(PupilFunction,allA,allB,allD,false);
            FieldMatrixDerivatives{itel,jtel,jz,lambdaderindex} = FieldMatrixDerivatives{itel,jtel,jz,lambdaderindex}...
              +((XImage_tmp-xemit)/lambda).*FieldMatrixDerivatives{itel,jtel,jz,1}...
              +((YImage_tmp-yemit)/lambda).*FieldMatrixDerivatives{itel,jtel,jz,2}...
              +((zemitrun-zemit)/lambda)*FieldMatrixDerivatives{itel,jtel,jz,zderindex};
        end

        % pupil functions for Zernike mode-derivative and FT to matrix elements
        switch parameters.fitmodel
          case 'aberrations'
            for jzer = 1:numzers
              jder = 3+jzer;
              PupilFunction = (2*pi*1i*normfac(jzer)*squeeze(allzernikes(jzer,:,:))/lambda).*PositionPhaseMask.*PupilMatrix{itel,jtel};
%               IntermediateImage = transpose(cztfunc(PupilFunction,Ay,By,Dy));
%               FieldMatrixDerivatives{itel,jtel,jz,jder} = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
              FieldMatrixDerivatives{itel,jtel,jz,jder} = cztfuncn(PupilFunction,allA,allB,allD,false);
            end
        end
      else
        for jder = 1:numders
          FieldMatrixDerivatives{itel,jtel,jz,jder} = zeros(Mx,My); % fill with dummy values
        end
      end
      
    end
  end
end

end
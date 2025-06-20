function [PSF,PSFderivatives] = get_psfs_derivatives(FieldMatrix,FieldMatrixDerivatives,parameters,compders)
% This function calculates the free or fixed dipole PSFs given the field
% matrix, the dipole orientation, and the pupil polarization, as well as
% the derivatives w.r.t. the xyz coordinates of the emitter and w.r.t. the
% emission wavelength lambda.
%
% copyright Sjoerd Stallinga, TU Delft, 2017
%
% parameters: emitter/absorber dipole orientation (characterized by angles
% pola and azim), detection/illumination polarization in objective lens
% back aperture (characterized by angles alpha and beta).
%
% flag compders indicates mode with/without computation PSFderivatives 

pola = parameters.pola;
azim = parameters.azim;
polarizationpupil = parameters.polarizationpupil;
alpha = parameters.alpha;
beta = parameters.beta;

dipor(1) = sin(pola)*cos(azim);
dipor(2) = sin(pola)*sin(azim);
dipor(3) = cos(pola);

polpupil(1) = cos(alpha)*exp(1i*beta);
polpupil(2) = sin(alpha)*exp(-1i*beta);

% find dimensions and number of derivatives from input
dims = size(FieldMatrixDerivatives);
if (length(dims)>3)
  Mz = dims(3);
  numders = dims(4);
  imdims = size(FieldMatrix{1,1,1});
else
  Mz = 1;
  numders = dims(3);
  imdims = size(FieldMatrix{1,1});
end
Mx = imdims(1);
My = imdims(2);

% calculation of free and fixed dipole PSF and the derivatives for the focal stack
PSF = zeros(Mx,My,Mz);
PSFderivatives = zeros(Mx,My,Mz,numders);
Ecder = cell(1,numders);
Exder = cell(1,numders);
Eyder = cell(1,numders);

for jz = 1:Mz
  
% calculation of free PSF and derivatives 
  if strcmp(parameters.dipoletype,'free')
    for jtel = 1:3
      if (polarizationpupil)
        Ec = polpupil(1)*FieldMatrix{1,jtel,jz}+polpupil(2)*FieldMatrix{2,jtel,jz};
        PSF(:,:,jz) = PSF(:,:,jz) + (1/3)*abs(Ec).^2;
      else
        for itel = 1:2
          PSF(:,:,jz) = PSF(:,:,jz) + (1/3)*abs(FieldMatrix{itel,jtel,jz}).^2;
        end
      end
      if compders
        for jder = 1:numders
          if (polarizationpupil)
            Ecder{jder} = polpupil(1)*FieldMatrixDerivatives{1,jtel,jz,jder}+polpupil(2)*FieldMatrixDerivatives{2,jtel,jz,jder};
            PSFderivatives(:,:,jz,jder) = PSFderivatives(:,:,jz,jder) + (2/3)*real(conj(Ec).*Ecder{jder});
          else
            for itel = 1:2
              PSFderivatives(:,:,jz,jder) = PSFderivatives(:,:,jz,jder) +...
                (2/3)*real(conj(FieldMatrix{itel,jtel,jz}).*FieldMatrixDerivatives{itel,jtel,jz,jder});
            end
          end
        end
      end
    end
  end
    
% calculation of fixed PSF and derivatives
  if strcmp(parameters.dipoletype,'fixed')
    Ex = dipor(1)*FieldMatrix{1,1,jz}+dipor(2)*FieldMatrix{1,2,jz}+dipor(3)*FieldMatrix{1,3,jz};
    Ey = dipor(1)*FieldMatrix{2,1,jz}+dipor(2)*FieldMatrix{2,2,jz}+dipor(3)*FieldMatrix{2,3,jz};
    if (polarizationpupil)
      Ec = polpupil(1)*Ex+polpupil(2)*Ey;
      PSF(:,:,jz) = abs(Ec).^2;
    else
      PSF(:,:,jz) = abs(Ex).^2+abs(Ey).^2;
    end
    if compders
      for jder = 1:numders
        Exder{jder} = dipor(1)*FieldMatrixDerivatives{1,1,jz,jder}+dipor(2)*FieldMatrixDerivatives{1,2,jz,jder}+dipor(3)*FieldMatrixDerivatives{1,3,jz,jder};
        Eyder{jder} = dipor(1)*FieldMatrixDerivatives{2,1,jz,jder}+dipor(2)*FieldMatrixDerivatives{2,2,jz,jder}+dipor(3)*FieldMatrixDerivatives{2,3,jz,jder};
        if (polarizationpupil)
          Ecder{jder} = polpupil(1)*Exder{jder}+polpupil(2)*Eyder{jder};
          PSFderivatives(:,:,jz,jder) = PSFderivatives(:,:,jz,jder) +...
            2*real(conj(Ec).*Ecder{jder});
         else
           PSFderivatives(:,:,jz,jder) = PSFderivatives(:,:,jz,jder) +...
            2*real(conj(Ex).*Exder{jder})+2*real(conj(Ey).*Eyder{jder});
        end
      end
    end
  end

% blurring due to non-zero pixel size, added 20180411
  PSF(:,:,jz) = do_pixel_blurring(PSF(:,:,jz),parameters);
  if compders
    for jder = 1:numders
      PSFderivatives(:,:,jz,jder) = do_pixel_blurring(PSFderivatives(:,:,jz,jder),parameters);
    end
  end
  
end

% 3D convolution of the PSFs and derivatives with a bead
if isfield(parameters,'bead')
  if parameters.bead == true
    bead = create3DBead(parameters);
    PSF = convn(bead,PSF,'same');
    [Mx,My,Mz] = size(PSF);
    if compders
      tempderivs = zeros(Mx,My,Mz,numders);
      for jder = 1:numders
        tempderivs(:,:,:,jder) = convn(bead,squeeze(PSFderivatives(:,:,:,jder)),'same');
      end
      PSFderivatives = tempderivs;
    end
  end
end
 

end


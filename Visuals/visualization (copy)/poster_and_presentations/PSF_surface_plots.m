% This script plots a 3D-surface plot of the central slice of a Gaussian
% and vectorial PSF.

% Generate Gaussian PSF.
x1 = -3:0.2:3;
x2 = -3:0.2:3;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];
gaussian_psf = mvnpdf(X);
gaussian_psf = reshape(gaussian_psf,length(x2),length(x1));

% Generate vectorial PSF without aberrations.

% Input
parameters = set_parameters_2D_psf_noaberrations;
N = 1;
bg = 0;

[wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(parameters);
[FieldMatrix,FieldMatrixDerivatives] = get_field_matrix_derivatives(parameters,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
[PSF,~] = get_psfs_derivatives(parameters,PupilMatrix,FieldMatrix,FieldMatrixDerivatives);
PSF = N*PSF+bg;

% Plot the PSFs.
%surf(x1,x2,gaussian_psf);
%surf(gaussian_psf);
surf(PSF);

addpath('C:\Users\idroste\Desktop\TUDelft\Code\diplib\share\DIPimage');
setenv('PATH',['C:\Users\idroste\Desktop\TUDelft\Code\diplib\bin',';',getenv('PATH')]);
dipshow(PSF,'lin');
colormap('hot');
diptruesize(400);
%colormap hot
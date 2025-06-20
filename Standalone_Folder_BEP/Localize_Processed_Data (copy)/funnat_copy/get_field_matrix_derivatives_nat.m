function [FieldMatrix,FieldMatrixDerivatives] = get_field_matrix_derivatives_nat(params,PupilMatrix,allzernikes,wavevector,wavevectorzimm,FOV_x,FOV_y)
% This function calculates the field matrix A_{jk}, which gives the j-th
% electric field component proportional to the k-th dipole vector
% component, as well as the derivatives of A_{jk} w.r.t. the xyz coordinates
% of the emitter and w.r.t. the emission wavelength lambda.
%
% parameters: NA, refractive indices of medium, wavelength (in nm),
% nominal emitter position (in nm) with z-position from
% cover slip-medium interface, spot footprint (in nm), axial range (in nm),
% sampling in pupil with (even), sampling in image plane (odd), sampling in
% axial direction

xemit = params.xemit;
yemit = params.yemit;
zemit = params.zemit;
Npupil = params.Npupil;
fitmodel = params.fitmodel;
numgammas = params.numgammas;
numgammas_fitted = params.numgammas_fitted;
aberrations = params.aberrations;
fitted_gamma_indices = params.fitted_gamma_indices;
lambda = params.lambda;

% calculate auxiliary vectors for chirpz
Ax = params.Axmt;
Bx = params.Bxmt;
Dx = params.Dxmt;
Ay = params.Aymt;
By = params.Bymt;
Dy = params.Dymt;

% calculation Zernike mode normalization
orders = params.aberrations(:,1:2);
zernike_normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));

% Number of derivatives
numders = numgammas_fitted;

if strcmp(params.excitation,'zstack')
    error('params.excitation = zstack is not implemented for global fiting.')
    
else %strcmp(params.excitation,'zstack')
    
    % phase contribution due to position of the emitter
    Wlateral = xemit*wavevector(:,:,1)+yemit*wavevector(:,:,2);
    if strcmp(params.ztype,'stage')
        Wpos = Wlateral+zemit*wavevectorzimm;
    elseif strcmp(params.ztype,'medium')
        Wpos = Wlateral+zemit*wavevector(:,:,3);
    end
    PositionPhaseMask = exp(-1i*Wpos);
    
    % Pupil function and field matrix
    PupilFunction = PositionPhaseMask.*PupilMatrix;
    IntermediateImage = cztfunc2D(PupilFunction,Ay,By,Dy,params);
    FieldMatrix = cztfunc2D(IntermediateImage,Ax,Bx,Dx,params);

    % pupil function for gamma-derivatives
    orders2D = params.orders2D;
    alllegendres2D = get_legendrepolynomials2D(orders2D,FOV_x,FOV_y);
    legendre_normfac = params.legendre_normfac;

    PupilFunctionDerivatives = zeros(Npupil,Npupil,numders,2,3);
    
    if strcmp(fitmodel,'xy-gamma')
        jv = 1;
    elseif strcmp(fitmodel,'xyz-gamma')
        jv = 6;
        %jv = 1;% In this way, you cannot exclude other gammas. But it works for excluding defocus for 3D fitting.
    end
    for jder = 1:numders
        jgam = fitted_gamma_indices(jder);
        RR = params.NATgammas{jgam};
        multiples = length(RR)-1;
        for jm = 1:multiples
            NATvec = RR{jm};
            nn = NATvec(3);
            mm = NATvec(4);
            weight = NATvec(5);
            jzer = intersect(find(aberrations(:,1)==nn),find(aberrations(:,2)==mm)); % find the right Zernike in the list
            derivative_jm = (2*pi*1i/lambda)*weight*legendre_normfac(jv)*alllegendres2D(:,:,jv)*zernike_normfac(jzer)*allzernikes(:,:,jzer).*PupilFunction;
            PupilFunctionDerivatives(:,:,jder,:,:) = PupilFunctionDerivatives(:,:,jder,:,:) + ...
            reshape(derivative_jm, [size(derivative_jm,1), size(derivative_jm,2), 1, size(derivative_jm,3), size(derivative_jm,4)]);
            jv = jv+1;
        end
    end

    IntermediateImage = cztfunc3D_nat(PupilFunctionDerivatives,Ay,By,Dy,params);
    FieldMatrixDerivatives = cztfunc3D_nat(IntermediateImage,Ax,Bx,Dx,params);
end
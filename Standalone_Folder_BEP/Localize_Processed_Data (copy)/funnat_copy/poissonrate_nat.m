function [mu,dmudtheta] = poissonrate_nat(params,theta,PupilMatrix,allzernikes,wavevector,wavevectorzimm,FOV_x,FOV_y)
% Returns the Poisson-rates for all pixels and all first order derivatives
% w.r.t. the parameters theta.

K = params.K;
Mx = params.Mx;
My = params.My;

fitmodel = params.fitmodel;
numgammas_fitted = params.numgammas_fitted;

switch fitmodel
    case {'xy-gamma'}
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
    case {'xyz-gamma'}
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
    otherwise
        error('fitmodel is not correct or not implemented.')
end

switch params.excitation
    case 'constant'
        P = 1;
    case 'zstack'
        P = ones(1,K);
end

% update pupil function
if contains(fitmodel,'aberrations')
    [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);
end

[FieldMatrix,FieldMatrixDerivatives] = get_field_matrix_derivatives_nat(params,PupilMatrix,allzernikes,wavevector,wavevectorzimm,FOV_x,FOV_y);
[PSF,PSFder] = get_psfs_derivatives(params,PupilMatrix,FieldMatrix,FieldMatrixDerivatives);

% get Poisson rate and derivatives
mu = zeros(Mx,My,K);
dmudtheta = zeros(Mx,My,K,numgammas_fitted);

if strcmp(params.excitation,'zstack')
    error('params.excitation = zstack is not implemented for global fiting.')
else
    for k = 1:K
        % PSF model
        mu(:,:,k) = Nph*P(k)*PSF+Nbg/K;
        dmudtheta(:,:,k,:) = Nph*P(k)*PSFder(:,:,:);
    end
end
function [CRLB_surface,Fisher,Fisher_inverse,CRLB_gammas,CRLB_Wsys,rcondstore] = get_fisher_crlb_global(mustore,dmudthetastore,Xn,Yn,outliers,theta_global,params)
% fprintf('\nFisher information: '); tic;
% This function calculates the Fisher-matrix and the Cramer-Rao Lower Bound
% for the parameters found.

% [Xn, Yn] = meshgrid for CRLB surface.

keps = 1e3*eps;
Ncfg = params.Ncfg;
numdirs = params.numgammas_fitted;
numgammas = params.numgammas;
gammas_fitted = params.gammas_fitted;
numgammas_fitted = params.numgammas_fitted;
aberrations = params.aberrations;

Fisher = zeros(numdirs,numdirs);
size_grid = size(Xn);
nr_of_zernike_modes = size(params.aberrations,1);
CRLB_surface = zeros(size_grid(1),size_grid(2),nr_of_zernike_modes);

% Calculate Fisher matrix
for jcfg = setdiff(1:Ncfg,outliers)

    mu = mustore(:,:,:,jcfg);
    dmudtheta = dmudthetastore(:,:,:,:,jcfg);

    % calculation Poisson rates
    mu = mu+params.readnoisevariance;
    mupos = double(mu>0).*mu + double(mu<0)*keps;
    weight = 1./mupos;

    % calculation Fisher matrix
    for ii = 1:numdirs
        for jj = ii:numdirs
            Fisher(ii,jj) = Fisher(ii,jj) + sum(weight.*dmudtheta(:,:,:,ii).*dmudtheta(:,:,:,jj),1:3);
            Fisher(jj,ii) = Fisher(ii,jj);
        end
    end
end

% regularization Fisher-matrix in order to circumvent possibility for
% inverting ill-conditioned matrix
% if (cond(Fisher)^-1>keps)
%     Fisher_inverse = inv(Fisher+keps*eye(size(Fisher)));
% end
% old version seems not correct. It only uses regularization when matrix is
% well conditioned

if (cond(Fisher)^-1>keps)
    Fisher_inverse_temp = inv(Fisher);
else
     Fisher_inverse_temp = inv(Fisher+keps*eye(size(Fisher)));
end
rcondstore = rcond(Fisher);

Fisher_inverse = zeros(numgammas,numgammas);
gammas_fitted_idx = find(gammas_fitted);
for i = 1:size(Fisher_inverse_temp,1)
    for j = 1:size(Fisher_inverse_temp,2)
        k = gammas_fitted_idx(i);
        m = gammas_fitted_idx(j);
        Fisher_inverse(k,m) = Fisher_inverse_temp(i,j);
    end
end

numders = numgammas;
daberrationdgamma = zeros(size_grid(1),size_grid(2),nr_of_zernike_modes,numgammas);

orders2D = params.orders2D;
alllegendres2D = get_legendrepolynomials2D(orders2D,Xn,Yn);
legendre_normfac = params.legendre_normfac;
jv = 1;
for jgam=1:numgammas
    RR = params.NATgammas{jgam};
    multiples = length(RR)-1;
    for jm = 1:multiples
        NATvec = RR{jm};
        nn = NATvec(3);
        mm = NATvec(4);
        weight = NATvec(5);
        jzer = intersect(find(aberrations(:,1)==nn),find(aberrations(:,2)==mm));
        daberrationdgamma(:,:,jzer,jgam) = weight*legendre_normfac(jv)*alllegendres2D(:,:,jv);
        jv = jv+1;
    end
end

% Calculate CRLB surface on meshgrid Xn,Yn 
for j = 1:numders
    for k = 1:numders
        for n=1:1:nr_of_zernike_modes
            CRLB_surface(:,:,n) = CRLB_surface(:,:,n) + Fisher_inverse(j,k).*daberrationdgamma(:,:,n,j).*daberrationdgamma(:,:,n,k);
        end
    end
end

CRLB_surface  = sqrt(CRLB_surface);

% Calculate CRLB of gammas
CRLB_gammas = sqrt(diag(Fisher_inverse));

% Calculate CRLB of system aberration value (Wsys)
CRLB_Wsys = 0;
multiplex_factors = params.multiplex_factors;
system_aberration_value = get_system_aberration_value(theta_global,params)*params.lambda/1e3; % in nm just like the other params.
for j = 1:numders
    for k = 1:numders
            CRLB_Wsys = CRLB_Wsys + (Fisher_inverse(j,k)*multiplex_factors(j)*multiplex_factors(k)*theta_global(j)*theta_global(k))/(system_aberration_value.^2);                                                                                                                                                                                                                                                                                                                                                                                                                                                          ;
    end
end
CRLB_Wsys = sqrt(CRLB_Wsys);

end

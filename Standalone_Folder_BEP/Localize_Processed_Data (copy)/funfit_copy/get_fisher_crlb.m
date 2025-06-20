function [CRLBstore,rcondstore,Fisher] = get_fisher_crlb(params,mustore,dmudthetastore)
% fprintf('\nFisher information: '); tic;
% This function calculates the Fisher-matrix and the Cramer-Rao Lower Bound
% for the parameters found.

keps = 1e3*eps;
Ncfg = params.Ncfg;

numdirs = params.numparams;
CRLBstore = zeros(numdirs,Ncfg);
rcondstore = zeros(1,Ncfg);
for jcfg = 1:Ncfg

    mu = mustore(:,:,:,jcfg);
    dmudtheta = dmudthetastore(:,:,:,:,jcfg);

    % calculation Poisson rates
    mu = mu+params.readnoisevariance;
    mupos = max(mu, keps);
    weight = 1./mupos;

    % calculation Fisher matrix
    Fisher = zeros(numdirs,numdirs);
    for ii = 1:numdirs
        for jj = ii:numdirs
            Fisher(ii,jj) = sum(weight.*dmudtheta(:,:,:,ii).*dmudtheta(:,:,:,jj),1:3);
            Fisher(jj,ii) = Fisher(ii,jj);
        end
    end

    % regularization Fisher-matrix in order to circumvent possibility for
    % inverting ill-conditioned matrix
    if (cond(Fisher)^-1>keps)
        CRLBstore(:,jcfg) = sqrt(diag(inv(Fisher+keps*eye(size(Fisher)))));
    end
    rcondstore(jcfg) = rcond(Fisher);
end

end

function [thetanew, flip_z, invertible] = thetaupdate_nat(thetaold,thetamax,thetamin,thetaretry,grad,Hessian,alambda,params)
% This function calculates the new parameters in each iteration step.

invertible = true;
numgammas_fitted = params.numgammas_fitted;
fitmodel = params.fitmodel;
flip_z = false;

% update of fit parameters via Levenberg-Marquardt

Bmat = Hessian+alambda*diag(diag(Hessian));

% Compute the inverse using "\". This uses gaussian elimination. If
% the matrix is not invertible, use a pseudoinverse. This requires
% less computations than always using the pseuodinverse.
if (abs(det(Bmat))>2*eps)
    dtheta = -Bmat\grad';
else
    invertible = false;
    dtheta = pinv(-Bmat)*transpose(grad);
end

thetanew = thetaold;
for i=1:numgammas_fitted
    j = params.fitted_gamma_indices(i);
    thetanew(j) = thetaold(j) + dtheta(i);
end
% 
if thetanew(5)<0
    thetanew = thetanew.*[-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1]';
    if strcmp(fitmodel,'xyz-gamma')
        flip_z = true;
    end
end

% enforce physical boundaries in parameter space.
for jj=1:length(thetaold)
    if ((thetanew(jj)>thetamax(jj))||(thetanew(jj)<thetamin(jj)))
        thetanew(jj) = thetaretry(jj);
    end
end

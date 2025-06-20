function [merit_total,merit_allspots,grad_total,Hessian_total,dmudtheta_allspots] = get_merit_nat(theta_local,theta_global,allspots,roixy,params)
%get_merit Calculates the total merit for a given choice of the local and
%global variables, for all spots.
% The total merit is the sum of the merit over all spots.

Ncfg = params.Ncfg;
Mx = params.Mx;
My = params.My;
K = params.K;
selected_indices_global = params.selected_indices_global;

varfit = params.varfit;
numgammas_fitted = params.numgammas_fitted;
dmudtheta_allspots = zeros(Mx,My,K,numgammas_fitted,Ncfg);

% Compute zernike coefficients
fov_coordinates = get_fov_coordinates(roixy,theta_local(1,:),theta_local(2,:),params);
fov_x = fov_coordinates(1,:);
fov_y = fov_coordinates(2,:);
zernikeCoefficients = get_zernike_coefficients(fov_x, fov_y, theta_global, params);

merit_allspots = zeros(1,Ncfg);
grad_allspots = zeros(numgammas_fitted,Ncfg);
Hessian_allspots = zeros(numgammas_fitted,numgammas_fitted,Ncfg);

% Evaluate merit
parfor jcfg=1:Ncfg
%for jcfg=1:Ncfg
    if ismember(jcfg,selected_indices_global)
        [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params,zernikeCoefficients(jcfg,:));
        [mu,dmudtheta] = poissonrate_nat(params,theta_local(:,jcfg),PupilMatrix,allzernikes,wavevector,wavevectorzimm,fov_x(jcfg),fov_y(jcfg));
        dmudtheta_allspots(:,:,:,:,jcfg) = dmudtheta;
    
        spot = allspots(:,:,:,jcfg);
        [merit,grad,Hessian] = likelihood_nat(params,spot,mu,dmudtheta,varfit);
        merit_allspots(jcfg) = merit;
        
        % summing outside of loop
        grad_allspots(:,jcfg) = grad;
        Hessian_allspots(:,:,jcfg) = Hessian;
    end
end

grad_total = sum(grad_allspots,2)';
Hessian_total = sum(Hessian_allspots,3);
merit_total = sum(merit_allspots,"all");

end
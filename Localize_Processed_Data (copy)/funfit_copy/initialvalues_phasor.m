function [thetainit] = initialvalues_phasor(allspots,roixy,theta_global,params)
% This function provides initial values for the fit parameters using the
% phasor method and a least squares fit.
%
% First, the initial position (x0,y0) is determined using the phasor method.
% After that, an initial estimate is made for the background b0 by taking 
% the median of the rim pixels and for the total photon count N0 by taking 
% the total photon of the measured spot minus b0.
% Then, this initial estimate for b0 and N0 is improved using by linear
% regression, i.e. a linear least squares fit to n = N*PSF+b.
% z0 is set to 0.
%
% Author: Sjoerd Stallinga, 2017-2021, Isabel Droste, 2022-2024, TU Delft
% Based on: Phasor based single-molecule localization microscopy in 3D
% (pSMLM-3D): An algorithm for MHz localization rates using standard CPUs
% Koen J. A. Martens, Arjen N. Bader, Sander Baas, Bernd Rieger, Johannes Hohlbein 
%
% Input
% allspots: all the measured spots in array of size (Mx, My, Mz, Ncfg).
% params: various parameters as specified in the set_parameters file.
%
% Output
% thetainit: The initial estimate of the fit parameters for all spots. This
% is an array of size (numparams,Ncfg).

K = params.K;
Mx = params.Mx;
My = params.My;
Ncfg = params.Ncfg;
fitmodel = params.fitmodel;
numparams = params.numparams;
pixelsize = params.pixelsize;
excitation = params.excitation;

% Get a meshgrid of physical coordinates in image plane.
[~,~,XX,YY] = get_coords(params);

% Allocate memory.
thetainit_x = zeros(1,Ncfg);
thetainit_y = zeros(1,Ncfg);
thetainit_z = zeros(1,Ncfg);
thetainit_N = zeros(1,Ncfg);
thetainit_b = zeros(1,Ncfg);
thetainit_aberrations = zeros(numel(params.aberrations(:,3)),Ncfg);

%%%%%%% new
% Calculate zernike coefficients to improve the initial values. 
%if strcmp(fitmodel,'xy-gamma') || strcmp(fitmodel,'xyz-gamma')
fov_coordinates = get_fov_coordinates(roixy,zeros(1,Ncfg),zeros(1,Ncfg),params);
zernikeCoefficients = get_zernike_coefficients(fov_coordinates(1,:), fov_coordinates(2,:), theta_global, params);
%else
%    zernikeCoefficients = repmat(params.aberrations(:,3),[1,Ncfg])';
%end
%%%%%%%%%%%%%

% Set up parallel pool
if isempty(gcp('nocreate'))
    number_of_parallel_workers = params.nr_of_parallel_workers;
    parpool('threads', number_of_parallel_workers);
end

% Compute the initial values for all spots.
parfor jcfg = 1:Ncfg
%parfor jcfg = 1:1000
%for jcfg = 1:Ncfg
    if rem(jcfg,round(Ncfg/10)) == 0
        fprintf('initial value %i\n',jcfg);
    end

    % Select spot
    fxyk = allspots(:,:,:,jcfg);
    if strcmp(excitation,'zstack')
        fxy = fxyk(:,:,round(K/2));
    else
        fxy = sum(fxyk,3); %why do we sum when it is not a z-stack?
    end

    % calculation of the phasor components for estimating the lateral
    % position, following Martens et al, JChemPhys 2018
    % this is more robust against background noise than centroid and
    % results in almost the right values
    qx1 = 2*pi/Mx/pixelsize;
    qy1 = 2*pi/My/pixelsize;
    sinfacx = sin(qx1*XX);
    cosfacx = cos(qx1*XX);
    sinfacy = sin(qy1*YY);
    cosfacy = cos(qy1*YY);
    phasor_sinx = sum(sum(sinfacx.*fxy));
    phasor_cosx = sum(sum(cosfacx.*fxy));
    phasor_siny = sum(sum(sinfacy.*fxy));
    phasor_cosy = sum(sum(cosfacy.*fxy));
    x0 = atan2(phasor_sinx,phasor_cosx)/qx1;
    y0 = atan2(phasor_siny,phasor_cosy)/qy1;

    % initial estimate of background bg from median value of the rim pixels,
    % and signal photon count Nph from total photon count in ROI
    rimpixels = zeros(2*Mx+2*My-4,1);
    rimpixels(1:Mx-1) = fxy(1:Mx-1,1);
    rimpixels(Mx:2*Mx-2) = fxy(2:Mx,My);
    rimpixels(2*Mx-1:2*Mx+My-3) = fxy(Mx,1:My-1);
    rimpixels(2*Mx+My-2:2*Mx+2*My-4) = fxy(1,2:My);
    bg0 = median(rimpixels);
    bg0 = max(bg0,1);
    Nph0 = sum(sum(fxy))-Mx*My*bg0;
    if (Nph0<0)
        Nph0 = sum(sum(fxy));
    end
    
    % z0 estimator
    if contains(fitmodel,'xyz') 
    aberrations = params.aberrations;
    jzer_horver = intersect(find(aberrations(:,1)==2),find(aberrations(:,2)==2)); % find the right Zernike in the list
    if isempty(jzer_horver)
      astighorver = 0.0;
    else
      astighorver = zernikeCoefficients(jcfg,jzer_horver)/params.lambda;
    end
    jzer_diag = intersect(find(aberrations(:,1)==2),find(aberrations(:,2)==-2)); % find the right Zernike in the list
    if isempty(jzer_diag)
      astigdiag = 0.0;
    else
      astigdiag = zernikeCoefficients(jcfg,jzer_diag)/params.lambda;
    end
    astigrms = sqrt(astighorver^2+astigdiag^2);
    astangle = atan2(astigdiag,astighorver)/2;
    % [astigrms_nm,idx] = max([zernikeCoefficients(jcfg,jzer_diag),zernikeCoefficients(jcfg,jzer_vert)]);
    % astigrms = astigrms_nm/params.lambda;
    % if idx==1
    %     astangle = 45*pi/180; % diagonal astigmatism
    % else
    %     astangle = 0; % horizontal/vertical astigmatism
    % end
    sinfacxa = sin(qx1*(XX-x0));
    sinfacya = sin(qy1*(YY-y0));
    cosfacxa = cos(qx1*(XX-x0));
    cosfacya = cos(qy1*(YY-y0));
    phasora_xy = sum(sum(sinfacxa.*sinfacya.*fxy));
    phasora_xx = 2*sum(sum(cosfacxa.*fxy));
    phasora_yy = 2*sum(sum(cosfacya.*fxy));
    deltA = cos(2*astangle)*(phasora_yy-phasora_xx)+sin(2*astangle)*2*phasora_xy;
    Ac = phasora_xx+phasora_yy;
    if params.NA<params.refmed
        dof = params.lambda/(params.refmed-sqrt(params.refmed^2-params.NA^2)); % depth of focus
    else
        dof = params.lambda/params.refmed;
    end
    scurvelength = (2*params.refmed*params.lambda/params.NA^2)*sqrt(8)*astigrms; % half the distance between focal lines, paraxial formula 
    gamfac = dof/scurvelength;
    if Ac^2>(1+gamfac^2)*deltA^2
        z0 = scurvelength*(1+gamfac^2)*deltA/(Ac+sqrt(Ac^2-(1+gamfac^2)*deltA^2));
    else
        z0 = scurvelength*(1+gamfac^2)*deltA/Ac;
    end
    else
        z0 = 0;
    end
    % Check that x0,y0,z0,are within limits
    if contains(fitmodel,'xyz') 
        thetatemp = [x0 y0 z0 Nph0 bg0];
    else 
        thetatemp = [x0 y0 Nph0 bg0];
    end
    [thetamin,thetamax] = thetalimits(params,thetatemp);
    
    
    for jj=1:length(thetatemp)-2
        % enforce physical boundaries in parameter space by 'thrust region reflective opteration'
        if thetatemp(jj) > thetamax(jj)
            thetatemp(jj) = 2*thetamax(jj) - thetatemp(jj);
        end
        if thetatemp(jj) < thetamin(jj)
            thetatemp(jj) = 2*thetamin(jj) - thetatemp(jj);
        end
        % double check to enforce physical boundaries in parameter space
        if (thetatemp(jj) > thetamax(jj) || thetatemp(jj) < thetamin(jj))
            thetatemp(jj) = (thetamin(jj)+thetamax(jj))/2;
        end
    end

    if contains(fitmodel,'xyz') 
        x0 = thetatemp(1);
        y0 = thetatemp(2);
        z0 = thetatemp(3);
        Nph0 = thetatemp(4);
        bg0 = thetatemp(5);
    else 
        x0 = thetatemp(1);
        y0 = thetatemp(2);
        Nph0 = thetatemp(3);
        bg0 = thetatemp(4);
    end

    % refinement of initial estimate of photon count and background by linear
    % regression, i.e. a linear least squares fit to n = Nph*PSF+bg
    % this procedure is entirely free from ad-hoc parameters and gives a bias
    % free initial estimate
    
    % first we compute the PSF for this initial position estimate (x0,y0)
    params_tmp = params;
    params_tmp.xemit = x0;
    params_tmp.yemit = y0;
    params_tmp.zemit = z0;
    params_tmp.K = 1;
    [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = ...
      get_pupil_matrix(params_tmp,zernikeCoefficients(jcfg,:)); %%%% new : extra input value
    [FieldMatrix,FieldMatrixDerivatives] = ...
      get_field_matrix_derivatives(params_tmp,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
    [PSF,~] = get_psfs_derivatives(params_tmp,PupilMatrix,FieldMatrix,FieldMatrixDerivatives);

    % use this for a refinement with weighted linear least squares
    % which is here mathematically equivalent to solving the MLE equations
    % d(log(L))/dNph=0 and d(log(L))/dbg=0 for Nph and bg
    % this gives an additional improvement in number of MLE iterations
    weight = 1./(Nph0*PSF+bg0);
    % weight = ones(size(PSF));
    Hav = mean(weight(:).*PSF(:));
    H2av = mean(weight(:).*PSF(:).^2);
    nav = mean(weight(:).*fxy(:));
    nHav = mean(weight(:).*fxy(:).*PSF(:));
    wav = mean(weight(:));
    detty = H2av*wav-Hav^2;
    Nph = (nHav*wav-nav*Hav)/detty;
    bg = (nav*H2av-nHav*Hav)/detty;

    % provision for negative or too large estimates
    if (bg < 1 || bg > 5*bg0)
        Nph = Nph+Mx*My*(bg-bg0);
        bg = bg0;
    end
    if (Nph < 0 || Nph > 2*sum(sum(fxy)))
        Nph = Nph0;
    end

    switch fitmodel
        case {'xy','xy-gamma'}
            thetainit_x(jcfg) = x0;
            thetainit_y(jcfg) = y0;
            thetainit_N(jcfg) = Nph;
            thetainit_b(jcfg) = bg;

        case {'xyz','xyz-gamma'}
            thetainit_x(jcfg) = x0;
            thetainit_y(jcfg) = y0;
            thetainit_z(jcfg) = z0;
            thetainit_N(jcfg) = Nph;
            thetainit_b(jcfg) = bg;       
        case 'xyz-aberrations'
            thetainit_x(jcfg) = x0;
            thetainit_y(jcfg) = y0;
            thetainit_z(jcfg) = z0;
            thetainit_N(jcfg) = Nph;
            thetainit_b(jcfg) = bg;
            thetainit_aberrations(:,jcfg) = params.aberrations(:,3);
        case 'xyz-tilt-aberrations'
            thetainit_x(jcfg) = x0;
            thetainit_y(jcfg) = y0;
            thetainit_z(jcfg) = z0;
            thetainit_N(jcfg) = Nph;
            thetainit_b(jcfg) = bg;
            thetainit_tiltx(jcfg) = params.tilt(1);
            thetainit_tilty(jcfg) = params.tilt(2);
            thetainit_aberrations(:,jcfg) = params.aberrations(:,3);
        otherwise
            error('fitmodel incorrect or not implemented for phasor initial values.')
    end

end % End for loop over spots

thetainit = zeros(numparams,Ncfg);

switch fitmodel
    case {'xy','xy-gamma'}
        thetainit(1,:) = thetainit_x;
        thetainit(2,:) = thetainit_y;
        thetainit(3,:) = thetainit_N;
        thetainit(4,:) = thetainit_b;
    case {'xyz','xyz-gamma'}
        thetainit(1,:) = thetainit_x;
        thetainit(2,:) = thetainit_y;
        thetainit(3,:) = thetainit_z;
        thetainit(4,:) = thetainit_N;
        thetainit(5,:) = thetainit_b;     
    case 'xyz-aberrations'
        thetainit(1,:) = thetainit_x;
        thetainit(2,:) = thetainit_y;
        thetainit(3,:) = thetainit_z;
        thetainit(4,:) = thetainit_N;
        thetainit(5,:) = thetainit_b;
        thetainit(6:end,:) = thetainit_aberrations;
    case 'xyz-tilt-aberrations'
        thetainit(1,:) = thetainit_x;
        thetainit(2,:) = thetainit_y;
        thetainit(3,:) = thetainit_z;
        thetainit(4,:) = thetainit_N;
        thetainit(5,:) = thetainit_b;
        thetainit(6,:) = thetainit_tiltx;
        thetainit(7,:) = thetainit_tilty;
        thetainit(8:end,:) = thetainit_aberrations;

    otherwise
        error('fitmodel incorrect or not implemented for phasor initial values.')
end


end
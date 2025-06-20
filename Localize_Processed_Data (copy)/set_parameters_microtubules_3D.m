function params = set_parameters_microtubules_3D
% This function sets all parameters for vectorial PSF calculations

% flags
params.FlagOTF = false;
params.Mpsfx = 129;
params.Mpsfy = 129;
params.Mpsfz = 71;
params.Notfx = 48;
params.Notfy = 48;
params.Notfz = 20;
params.OTFGridSizeX = 10;
params.OTFGridSizeY = 10;

% c++ fitting
params.cpp = false;%true;
params.cpp_fitmode = "gpu-lowaccuracy";% "gpu-lowaccuracy"; %"cpu", "gpu" or "gpu-lowaccuracy" double quotes needed. 
% to do: rename gpu-lowaccuracy->gpu and gpu->gpu-doubleprecision
params.initialization = 'phasor'; % to do: in matlab there is no choice 

% parameters: NA, refractive indices of medium, cover slip, immersion fluid,
% nominal value of immersion fluid refractive index matching objective lens
% design, nominal free working distance (in nm), distance image plane from
% cover slip (in nm), wavelength (in nm), emitter position (in nm) with
% z-position from image plane, spot footprint (in nm), axial range min/max
% (in nm), flag for axial range by z-position in medium or by z-stage
% position, sampling in pupil with, sampling in image plane, sampling in
% axial direction.
%

params.datatype = 'tif'; %'mat' 'tif', 'h5' % datatype of the raw data
params.simulated_data = false;
params.show_plots = true;

% fitting parameters
params.nr_of_parallel_workers = feature('numcores');
params.tollim = 5e-10;
params.tollim_total = 1e-6; % Stopping parameter of total aberration fitting loop (mlambda)
params.varfit = 0;

params.select_spots = true; % Select subset of Ncfg_total for finding field dependent aberrations.
if params.select_spots
    params.Ncfg_global = 1000; % selection out of Ncfg for global update
    params.Ncfg = 5*params.Ncfg_global; % all spots for global update
else
    params.Ncfg_global = params.Ncfg_total; % selection out of Ncfg for global update
    params.Ncfg = params.Ncfg_total; % all spots for global update
end

% Segmentation parameters
% You can check in the figure if the threshold is OK
params.segmentation_ROI_thr = 15.5;
params.nr_of_cutoff_frames = 0;

params.fit_aberrations = true;
params.multiple_initial_gammas = false;
params.min_total_iterations = 5;
params.max_total_iterations = 20;
params.min_local_iterations = 1;
params.max_local_iterations = 20;
params.min_global_iterations = 1;
params.max_global_iterations = 20;

params.perform_final_iterations = true;
params.max_final_iterations = 10; % Switch small number of local and global iteration to converge spherical aberration
params.per_final_iteration = 2; % Number of iterations per final iteration

params.alambda0_local = 1; % initial value of alambda (Levenberg Marquardt parameter)
params.alambda0_global = 1e3;
params.alambdafac_local = 10; % multiplication factor (Levenberg Marquardt parameter)
params.alambdafac_global = 10;

% camera offset and gain
params.offset_image = false;
params.gain_image = false;
params.offset = 100; % [AUD]
params.gain = 2.12; %[AUD/e-]

% global image parameters
params.imgSizeX = 1024;
params.imgSizeY = 1024;

% PSF/optical parameters
params.NA = 1.35;
params.refmed = 1.38; % n_sample
params.refcov = 1.52; % n_coverslips
params.refimm = 1.4; % n_immersion fluid; Olympus Type F: n = 1.518
params.refimmnom = params.refcov; % n of norminal immersion medium of objective. params.refimmnom = params.refcov; %params.refimmnom = 1; % for air objective
params.fwd = 0.2e6; % Working distance in nm
params.depth = 0; % position of sample along z. [nm]
params.zrange = [-100*8,100*8]; % z-stack range in nm
params.zspread = [-1000,1000]; % min and max z-value
params.ztype =  'stage'; % 'medium'
params.lambda = 680; % emission central wavelength
params.lambdacentral = 680;
params.lambdaspread = [680,680];
params.xemit = 0.0; % position of emitter in ROI
params.yemit = 0.0;
params.zemit = 0.0;
params.Npupil = 48; % = params.Mx+1
params.pixelsize = 95; % nm/px of sample image
params.samplingdistance = params.pixelsize;
params.Mx = 17;% segmentation size. If already segmented, can set to the size of cropped image
params.My = params.Mx;
params.Mz = 1;% total frame of z-stack. params.Mz=params.k  % sampling in axial direction

params.xrange = params.pixelsize*params.Mx/2;
params.yrange = params.pixelsize*params.My/2;

params.debugmode = 0;
params.flg_parallel = 0;
params.flg_nat = 0;

% SAF check
if and(params.NA>params.refmed, params.depth<4*params.lambda)
    [zvals,~] = set_saffocus(params);
else
    [zvals,~] = get_rimismatchpars(params);
end
params.zvals = zvals;

% sanity check on position emitter w.r.t. cover slip
if strcmp(params.ztype,'stage')
    zcheck = params.depth+params.zemit;
end
if strcmp(params.ztype,'medium')
    zmin = params.zrange(1);
    zcheck = zmin+params.depth+params.zemit;
end
if (zcheck<0)
    fprintf('Warning! Emitter must be above the cover slip:\n')
    fprintf('Adjust parameter settings for physical results.\n')
end

% sanity check on refractive index values
if (params.NA>params.refimm)
    fprintf('Warning! Refractive index immersion medium cannot be smaller than NA.\n')
end
if (params.NA>params.refcov)
    fprintf('Warning! Refractive index cover slip cannot be smaller than NA.\n')
end

% parameters needed for fixed dipole PSF only: emitter/absorber dipole
% orientation (characterized by angles pola and azim)
params.dipoletype = 'free';
params.pola = 90.0*pi/180;
params.azim = 0.0*pi/180;

% diffusion coefficient
welldepth = eps;
g2 = (3+welldepth^2-3*welldepth*coth(welldepth))/welldepth^2;
params.welldepth = welldepth;
params.g2 = g2;

params.fitmodel = 'xyz-gamma';

if strcmp(params.fitmodel,'xy-gamma') || strcmp(params.fitmodel,'xyz-gamma')

    params.wrms_level = 0.036; % not used (you can change it in generate_NAT_data)
    params.diffraction_limited_fraction = 1.0;

    % aberration order, must be j=0 (only defocus/field curvature),
    % or j=1 (first order aberration theory), we still need to work out the
    % case j=2 (second order aberration theory)
    j = 1; 
    
    % Set up array with Zernike aberration coefficients
    % total number of Zernike's up to order j, excluding tip/tilt, derived
    % using Mathematica: Nzer = sum[2p+1,{p,0,J}]-2
    Nzer = (j+3)*(j+1)-2;
    params.aberrations = zeros(Nzer,3);
    jzer = 1;
    for n = 2:2*(j+1)
      mmin = max([-n,n-2*(j+1)]);
      mmax = min([n,2*(j+1)-n]);
      for m = mmin:2:mmax
        params.aberrations(jzer,1) = n;
        params.aberrations(jzer,2) = m;
        params.aberrations(jzer,3) = 0.0; % dummy value
        jzer = jzer+1;
      end
    end

    % Set up cell array with Legendre NAT coefficients, I have not found a
    % generic formulation for every aberration order, hence we use a switch
    % statement here, later we should evaluate the equations for 2nd order
    % aberration theory
    %
    % NATgammas is an Ngam x 1 cell array, as there are Ngam different
    % Legendre NAT coefficients. For each entry NATgammas{jgam} we obtain a
    % nested cell array RR of size (M+1) x 1. The last value is the value of the
    % Legendre NAT coefficient. The first M elements are a set of row vectors
    % of the form [px py n m w]. Here (n,m) are the radial and azimuthal
    % Zernike indices and (px,py) are the Legendre indices. The weight w is
    % introduced to properly take into account the NAT symmetries between the
    % different NAT coefficients. We call M the multiplicity of the Legendre
    % NAT coefficient. The Legendre NAT coefficients are the RMS values, where
    % we choose the convention that the system RMS aberration is given by
    % <W_rms^2> = sum( M_j*gamma_j^2, {j=1,jmax} ) with M_jthe multiplex
    % factor. In the example below all gamma's are set to zero.
    
    
    switch j
      case 0
        params.NATgammas = {{[0 0 2 0 1.0],0.0}}; % trivial case only included pro forma
      case 1
        params.NATgammas = {{[0 0 2 0 1.0],0.0},...
                     {[1 0 2 0 1.0],0.0},...
                     {[0 1 2 0 1.0],0.0},...
                     {[2 0 2 0, 1.0],[0 2 2 0, 1.0],0.0},...
                     {[0 0 2 2 1.0],0.0},...
                     {[0 0 2 -2 1.0],0.0},...
                     {[1 0 2 2 1.0],[0 1 2 -2 1.0],0.0},...
                     {[0 1 2 2 -1.0],[1 0 2 -2 1.0],0.0},...
                     {[2 0 2 2 1/sqrt(5)],[0 2 2 2 -1/sqrt(5)],[1 1 2 -2 1.0],0.0},...
                     {[0 0 3 1 1.0],0.0},...
                     {[0 0 3 -1 1.0],0.0},...
                     {[1 0 3 1 1.0],[0 1 3 -1 1.0],0.0},...
                     {[0 0 4 0 1.0],0.0}};
    end
    params.numgammas = length(params.NATgammas);
    % Indicate which gammas are fitted (true) and which are fixed
    % (false). Don't change it, does not work.
    params.gammas_fitted = [ ...
    true; %1 defocus
    true; %2 defocus
    true; %3 defocus
    true; %4 defocus
    true; %5 astigmatism
    true; %6 astigmatism
    true; %7 astigmatism
    true; %8 astigmatism
    true; %9 astigmatism
    true; %10 coma
    true; %11 coma
    true; %13 coma
    true; %14 spherical aberration
    ];
    if strcmp(params.fitmodel,'xyz-gamma')
        % Indicate which gammas are fitted (true) and which are fixed
        % (false).
        params.gammas_fitted(1:4) = false;
    end
    params.numgammas_fitted = sum(params.gammas_fitted == true);
    params.fitted_gamma_indices = find(params.gammas_fitted);

    params.system_aberration_level = 40; % in mlambda
    params.aberration_level_init = 1;

    params.gammas_init = zeros(params.numgammas,1);
    params.gammas_init(6) = 150*params.lambda*1e-3;

    % compute multiplex factors and numNATparvecs
    params.multiplex_factors = zeros(params.numgammas,1); 
    params.numNATparvecs = 0;
    for jgam = 1:params.numgammas
      RR = params.NATgammas{jgam};
      multiples = length(RR)-1;
      for jm = 1:multiples
        NATvec = RR{jm};
        params.multiplex_factors(jgam) = params.multiplex_factors(jgam) + NATvec(5)^2;
      end
      params.numNATparvecs = params.numNATparvecs + multiples;
    end

    % List legendre orders
    params.orders2D = zeros(params.numNATparvecs,2);
    jv = 1; % counter to step through list of all different NATvec terms
    for jgam = 1:params.numgammas
      RR = params.NATgammas{jgam};
      multiples = length(RR)-1;
      for jm = 1:multiples
        NATvec = RR{jm};
        params.orders2D(jv,1) = NATvec(1);
        params.orders2D(jv,2) = NATvec(2);
        jv = jv+1;
      end
    end

    params.legendre_normfac = sqrt((1+2*params.orders2D(:,1)).*(1+2*params.orders2D(:,2)));

end

% DOE/SLM
params.doetype = 'none';
% params.doetype = 'vortex';
% params.ringradius = 1;
% params.doelevels = 64;
% params.zonefunction = params.aberrations;
% params.doephasedepth = 593;
% params.doevortexflip = 1;

% Bead parameters for convolution with PSF and derivatives, beaddiameter in nm
params.bead = false;
params.beaddiameter = 100;
% check on meaningfullness of bead convolution
if params.beaddiameter<params.pixelsize
    params.bead = false;
end

% Fit model parameters: signal photon count, background photons/pixel, read
% noise variance for sCMOS camera's, fit model, output labels depending on
% fit modelnn
params.readnoisestd = 0;
params.readnoisevariance = 0;

% model parameters
params.alpha = 0;
params.beta = 0;
params.K = params.Mz;
params.m = 1;
params.excitation = 'constant';%'zstack';

if strcmp(params.fitmodel,'xy')
    params.numparams = 4;
elseif strcmp(params.fitmodel,'xyz')
    params.numparams = 5;
elseif strcmp(params.fitmodel,'xy-gamma')
    params.numparams = 4;
elseif strcmp(params.fitmodel,'xyz-gamma')
    params.numparams = 5;
elseif strcmp(params.fitmodel,'xyz-aberrations')
    params.numparams = 5+numel(params.aberrations(:,3));
elseif strcmp(params.fitmodel,'xy-azim')
    params.numparams = 5;
elseif strcmp(params.fitmodel,'xy-azim-pola')
    params.numparams = 6;
elseif strcmp(params.fitmodel,'xyz-azim-pola')
    params.numparams = 7;
elseif strcmp(params.fitmodel,'xy-azim-diffusion')
    params.numparams = 6;
elseif strcmp(params.fitmodel,'xy-azim-pola-diffusion')
    params.numparams = 7;
elseif strcmp(params.fitmodel,'xyz-azim-pola-diffusion')
    params.numparams = 8;
elseif strcmp(params.fitmodel,'xyz-azim-pola-aberrations')
    params.numparams = 7+numel(params.aberrations(:,3));
elseif strcmp(params.fitmodel,'xyz-azim-pola-diffusion-aberrations')
    params.numparams = 8+numel(params.aberrations(:,3));
end

% % calculate auxiliary vectors for chirpz
PupilSize = params.NA/params.lambda;
[Ax,Bx,Dx] = prechirpz(PupilSize,params.xrange,params.Npupil,params.Mx);
[Ay,By,Dy] = prechirpz(PupilSize,params.yrange,params.Npupil,params.My);

params.Axmt = repmat(Ax,params.Mx,1);
params.Bxmt = repmat(Bx,params.Mx,1);
params.Dxmt = repmat(Dx,params.Mx,1);
params.Aymt = repmat(Ay,params.Npupil,1);
params.Bymt = repmat(By,params.Npupil,1);
params.Dymt = repmat(Dy,params.Npupil,1);

params.cztN = params.Npupil;
params.cztM = params.Mx;
params.cztL = params.Npupil+params.Mx-1;

% show intermediate results for monitoring code
params.debugmode = 0;
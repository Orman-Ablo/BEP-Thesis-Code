function [thetastore,mustore,dmudthetastore,meritstore,numiters, ...
    iteration_path,alambdastore,time] = ...
localization(allspots,theta0,params)
% This function finds the parameters for a single 2D-image assuming
% vectorial PSF models.
%
fprintf('\nStart fitting %i instances:\n',params.Ncfg);
total_time = tic;

% parameter settings
Ncfg = params.Ncfg;
varfit = params.varfit;
tollim = params.tollim;
Nitermax = params.Nitermax;
numparams = params.numparams;
nr_of_parallel_workers = params.nr_of_parallel_workers;

% pre-allocation
%numiters = zeros(1,Ncfg);
numiters = -1*ones(1,Ncfg);
mustore = zeros(params.Mx,params.My,params.K,Ncfg);
dmudthetastore = zeros(params.Mx,params.My,params.K,numparams,Ncfg);
thetastore = zeros(numparams,Ncfg,Nitermax+1);
meritstore = zeros(Ncfg,Nitermax+1);
time_per_spot = zeros(Ncfg,1); % time that it takes to localize each iteration.
time_per_iteration = zeros(Ncfg,Nitermax); % time that each iteration takes.


% Save which path is taken in levenberg-marquardt method in each iteration
% 01 : dmerit < 0 (No improvement, update not accepted, alambda increased).
%                  Normal inverse used (A\b)
% 02 : dmerit >= 0 (Improvement, update accepted, alambda decreased)
%                  Normal inverse used (A\b)
% 11 : same as 01 but pseudoinverse was used
% 12 : same as 02 but pseudoinverse was used
iteration_path = zeros(Ncfg,Nitermax);
alambdastore = zeros(Ncfg,Nitermax+1);

flg_nat = params.flg_nat;
if flg_nat
    natPredictions = params.natPredictions;
else
    natPredictions = zeros(Ncfg,1);
end

% setup parallel loop
if params.flg_parallel
    p = gcp('nocreate');
    if isempty(p)
        parpool(nr_of_parallel_workers);
    end
    parforArg = Inf;
else
    parforArg = 0;
end

parfor (jcfg = 1:Ncfg, parforArg)
%for jcfg = 1:Ncfg
    fprintf('\nStart fitting instance %i\n',jcfg);
    time_per_spot_counter = tic; 

    if flg_nat
        [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params,natPredictions(jcfg,:));
    else
        [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);
    end

    % pre-allocate
    thetatemp = zeros(numparams,Nitermax+1);
    merittemp = zeros(1,Nitermax+1);
    iteration_path_temp = zeros(Nitermax,1);
    alambdatemp = zeros(Nitermax+1,1);
    time_per_iteration_temp = zeros(Nitermax,1);
    
    % initial values and max/min values
    theta = theta0(:,jcfg);
    spots = allspots(:,:,:,jcfg);
    
    [thetamin,thetamax] = thetalimits(params,theta);
    
    thetaretry = (thetamax+thetamin)/2;
    
    [mu,dmudtheta] = poissonrate(params,theta,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
    [merit,grad,Hessian] = likelihood(params,spots,mu,dmudtheta,varfit);

    alambda = params.alambda0; 
    alambdafac = params.alambdafac;
    
    thetatemp(:,1) = theta;
    merittemp(1) = merit;
    alambdatemp(1,:) = alambda;
    
    % start iteration loop
    iiter = 1;
    monitor = 2*tollim;
    
    %while ((iiter<=Nitermax) && (monitor>tollim))
    while iiter<=Nitermax
        time_per_iteration_counter = tic;
         % save iteration in which convercenge would be reached.
        if ((monitor<=tollim) && (numiters(jcfg)<0))
            numiters(jcfg) = iiter -1;
        end
        if mod(iiter,1) == 0
            fprintf('\nBead %i: iteration %i\n',jcfg, iiter);
        end
 
        % update parameters
        [thetatry, invertible] = thetaupdate(theta,thetamax,thetamin,thetaretry,grad,Hessian,alambda,params);
          
        % calculate update merit function
        [mu,dmudtheta] = poissonrate(params,thetatry,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
        [merittry,gradtry,Hessiantry] = likelihood(params,spots,mu,dmudtheta,varfit);
        dmerit = merittry-merit;
            
        % modify Levenberg-Marquardt parameter
        if (dmerit<0)
            % If there if no improvement, do not accept the update,
            % and increase alambda.
            if invertible % Save that we are in this case.
                iteration_path_temp(iiter,:) = 01;
            else
                iteration_path_temp(iiter,:) = 11;
            end

            alambda = alambdafac*alambda;
            alambdatemp(iiter+1,:) = alambda;
        else
            % If there is an improvement, accept the update and
            % decrease alambda.
            if invertible % Save that we are in this case.
                iteration_path_temp(iiter,:) = 02;
            else
                iteration_path_temp(iiter,:) = 12;
            end

            alambda = alambda/alambdafac;
            theta = thetatry;
            merit = merittry;
            grad = gradtry;
            Hessian = Hessiantry;
            monitor = abs(dmerit/merit);
            thetaretry = theta;
            iteration_path_temp(iiter,:) = 2;
            alambdatemp(iiter+1,:) = alambda;
        end

        % store values and update counter
        thetatemp(:,iiter+1) = theta;
        merittemp(iiter+1) = merit;
          
        time_per_iteration_temp(iiter) = toc(time_per_iteration_counter);
        iiter = iiter+1; % update counter 
    end % end of while iiter<=Nitermax
    
    % store values
    %numiters(jcfg) = iiter-1;
    
    for jiter = iiter+1:Nitermax+1
        merittemp(jiter) = merit;
        thetatemp(:,jiter) = theta;
    end
    
    mustore(:,:,:,jcfg) = mu;
    dmudthetastore(:,:,:,:,jcfg) = dmudtheta;
    thetastore(:,jcfg,:) = thetatemp;
    meritstore(jcfg,:) = merittemp;

    iteration_path(jcfg,:) = iteration_path_temp;
    alambdastore(jcfg,:) = alambdatemp;
    
    time_per_spot(jcfg) = toc(time_per_spot_counter);
    time_per_iteration(jcfg,:) = time_per_iteration_temp;
end

% add offset to merit function
meritoffset = meritoffsetcalc(allspots,params.varfit);
for jiter = 1:Nitermax+1
    meritstore(:,jiter) = meritstore(:,jiter)+meritoffset;
end

% print run time
%fprintf(['\nMLE fit routine (spot/second): ' num2str(toc,3) 's (' num2str(params.Ncfg/toc,5) ')\n'])
time.total_time = toc(total_time);
time.time_per_spot = time_per_spot;
time.time_per_iteration = time_per_iteration;

end

%%% tmp fun
% calculate theta estimation limits (max/min values)
function [thetamin,thetamax] = thetalimits(params,theta)

fitmodel = params.fitmodel;
zernikecoefsmax = 0.25*params.lambda*ones(1,size(params.aberrations,1));

roisizex = params.Mx*params.pixelsize;
roisizey = params.My*params.pixelsize;
xmin = -roisizex/2;
xmax = roisizex/2;
ymin = -roisizey/2;
ymax = roisizey/2;
zmin = params.zspread(1);
zmax = params.zspread(2);
azimmax = 2*pi;
polamax = pi;

if strcmp(fitmodel,'xy')
    thetamin = [xmin,ymin,theta(3)/10,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3)];
elseif strcmp(fitmodel,'xy-azim')
    thetamin = [xmin,ymin,theta(3)/10,0,0];
    thetamax = [xmax,ymax,roisizey/2,2*theta(3),theta(3),azimmax];
elseif strcmp(fitmodel,'xy-azim-pola')
    thetamin = [xmin,ymin,theta(3)/10,0,0,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3),azimmax,polamax];
elseif strcmp(fitmodel,'xy-azim-pola-diffusion')
    thetamin = [xmin,ymin,theta(3)/10,0,0,-0,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3),azimmax,polamax,1];
end

if strcmp(fitmodel,'xyz')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4)];
elseif strcmp(fitmodel,'xyz-azim-pola')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax];
elseif strcmp(fitmodel,'xy-azim-diffusion')
    thetamin = [xmin,ymin,theta(3)/10,0,0,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3),azimmax,1];
elseif strcmp(fitmodel,'xyz-azim-pola-diffusion')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0,0];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax,1];
elseif strcmp(fitmodel,'xyz-aberrations')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,-zernikecoefsmax];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),zernikecoefsmax];
elseif strcmp(fitmodel,'xyz-azim-pola-aberrations')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0,-zernikecoefsmax];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax,zernikecoefsmax];
elseif strcmp(fitmodel,'xyz-azim-pola-diffusion-aberrations')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0,0,-zernikecoefsmax];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax,1,zernikecoefsmax];
end

end

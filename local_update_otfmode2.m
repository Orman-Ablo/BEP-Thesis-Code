function [theta_local,thetastore_local,localizations,localizations_with_outliers,meritstore_local,alambda_local,mu_allspots,dmudtheta_allspots,Niters_local,outliers] = ...
         local_update_otfmode2(theta_local,thetastore_local,meritstore_local,alambda_local,allspots,roixy,iiter_total,Niters_local,flip_z_net,framelist,ID,OTF_3d,intensity,params)
% This function executes an update of the local variables, consisting of
% multiple local iterations.
%
% Autor: Isabel Droste, TU Delft, 2022, 2023


% When symmetric aberrations have been flipped an uneven nr. of times in
% previous global update, z-coordinates need to be flipped
if flip_z_net
    theta_local(3,:) = -1*theta_local(3,:);
end

% Settings
Ncfg = params.Ncfg;
numparams = params.numparams;
tollim = params.tollim;
varfit = params.varfit;
max_local_iterations = params.max_local_iterations;
min_local_iterations = params.min_local_iterations;
fitmodel = params.fitmodel;
Npatch = params.Npatch;
Xpatch =params.Xpatch;
Ypatch =params.Ypatch;

% Allocation and initialization
mu_allspots = zeros(params.Mx,params.My,params.K,Ncfg);
dmudtheta_allspots = zeros(params.Mx,params.My,params.K,numparams,Ncfg);
Niters_temp = zeros(1,Ncfg);
thetastore_local_temp = zeros(numparams,Ncfg);
meritstore_local_temp = zeros(Ncfg,1);

% Compute aberrations
%if strcmp(fitmodel,'xy-gamma') || strcmp(fitmodel,'xyz-gamma')


%else
%    zernikeCoefficients = repmat(params.aberrations(:,3),[1,Ncfg])';
%end
    
%% Local update loop
fprintf('\nStart local update %i\n',iiter_total);


if params.check_test
    PSF_store = NaN(params.Mx,params.My,max_local_iterations+1,Ncfg);
    pos_store = NaN(max_local_iterations+1,3,Ncfg);
    grad_store = NaN(Ncfg,max_local_iterations+1,numparams);

end





%parfor jcfg = 1:Ncfg
%parfor jcfg = 1:1000
parfor jcfg = 1:Ncfg

    if rem(jcfg,round(Ncfg/10)) == 0
        fprintf('fitting spot %i\n',jcfg);
    end
    
    % inital values
    theta = theta_local(:,jcfg);
    spot = allspots(:,:,:,jcfg);
    [thetamin,thetamax] = thetalimits(params,theta);
    thetaretry = (thetamax+thetamin)/2;
    
    % TODO: OTF method.
    %
    % note: there is an overlap with the last iteration of the global update,
    % both compute the PSF, this can be combined.
    
    
    % if params.OTF
    %     params_local = params;
    %     params_local.xemit = 0;
    %     params_local.yemit = 0;
    %     % Calculate center PSF
    %     [FieldMatrix,FieldMatrixDerivatives] = get_field_matrix_derivatives(params_local,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
    %     [PSF,~] = get_psfs_derivatives(params_local,PupilMatrix,FieldMatrix,FieldMatrixDerivatives);
    % 
    %     % Calculate OTF from PSF
    %     OTF = get_otf(PSF,params);
    % end
    %params.OTF = false;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% xy-limits
    roix = roixy(1,jcfg);
    roiy = roixy(2,jcfg);
    x_loc = theta(1)+roix*params.pixelsize;% keep the row constant, but use a variable for the columns
    y_loc = theta(2)+roiy*params.pixelsize;
    % Localizations (xy) = thetaUpdate (xy)+roixy*pixelsize


    %% for each spot, the OTF is found ONCE: This means that as long as the right localization file was used, NO points 'should' fall outside the bound
    x_patch = NaN;
    y_patch = NaN;
    for ix = 1:Npatch
        if (x_loc>=Xpatch(ix))&&(x_loc<=Xpatch(ix+1))
            x_patch = ix;
            break
        end
    end
    
    for iy = 1:Npatch
        if (y_loc>=Ypatch(iy))&&(y_loc<=Ypatch(iy+1))
            y_patch = iy;
            break
        end
    end
    if isnan(x_patch) || isnan(y_patch)
        error('Patch index not found! x_loc=%.2f, y_loc=%.2f', x_loc, y_loc);
    end
    % selecting the right OTF for the position
    OTF = OTF_3d(:,:,:,x_patch,y_patch);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mu,dmudtheta] = poissonrate_otfmode2(params,theta,OTF,intensity); % modified Poissonrate file, to work with OTF
    [merit,grad,Hessian] = likelihood(params,spot,mu,dmudtheta,varfit);
    meritprev = merit;

    if params.check_test
        PSF_store(:,:,1,jcfg) = mu;
        pos_store(1,:,jcfg) = theta(1:3);
        grad_store(jcfg,1,:) = grad;

    end

    
    % start iteration loop
    iiter = 1;
    monitor = 2*tollim;
    alambda = alambda_local(jcfg);
    alambdafac = params.alambdafac_local;
    
    while ((iiter<=max_local_iterations) && ((iiter<=min_local_iterations) || (monitor>tollim)))

        %fprintf('spot %i iter %i\n',jcfg,iiter);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        % update parameters
        [thetatry,~] = thetaupdate(theta,thetamax,thetamin,thetaretry,grad,Hessian,alambda,params);

        % calculate update merit function
        [mutry,dmudthetatry] = poissonrate_otfmode2(params,thetatry,OTF,intensity);
        [merittry,gradtry,Hessiantry] = likelihood(params,spot,mutry,dmudthetatry,varfit);
        dmerit = merittry-merit;
        
        % modify Levenberg-Marquardt parameter
        if (dmerit<0)
            alambda = alambdafac*alambda;
        else
            alambda = alambda/alambdafac;
            theta = thetatry;
            mu = mutry;
            dmudtheta = dmudthetatry;
            merit = merittry;
            grad = gradtry;
            Hessian = Hessiantry;
            dmerit = merit-meritprev;
            monitor = abs(dmerit/merit);
            meritprev = merit;
            thetaretry = theta;
        end

        % if params.check_test
        %     PSF_store(:,:,iiter+1,jcfg) = mu;
        %     pos_store(iiter+1,:,jcfg) = theta(1:3);
        %     grad_store(jcfg,iiter+1,:) = grad;
        % end 
    
        % Store values and update counter
        thetastore_local_temp(:,jcfg) = theta;
        meritstore_local_temp(jcfg) = merit;

        iiter = iiter+1;
        
    end % end while iiter<max_local_iterations && monitor>tollim

    Niters_temp(jcfg) = iiter-1;
    if monitor>=tollim
        Niters_temp(jcfg) = iiter; 
    end

    theta_local(:,jcfg) = theta;
    alambda_local(jcfg) = alambda;     
    mu_allspots(:,:,:,jcfg) = mu;
    dmudtheta_allspots(:,:,:,:,jcfg) = dmudtheta;
    
end % end for jcfg = 1:Ncfg

Niters_local(:,iiter_total) = Niters_temp';
thetastore_local(:,:,end+1) = thetastore_local_temp;
meritstore_local(:,end+1) = meritstore_local_temp;
% Outliers are computed again each update, so the previous outliers are not
% remembered.
outliers = find(get_outliers(theta_local,meritstore_local(:,end),Niters_temp,mu_allspots,allspots,params));

% Calculate CRLB
CRLB_local = get_fisher_crlb(params,mu_allspots,dmudtheta_allspots);

% Create localization list in physical FOV coordinates. Format:
% x(nm) y(nm) (z(nm)) framelist CRLBx CRLBy CRLBz Nph bg
Ncfg_total = size(theta_local,2);
localizations_with_outliers = zeros(Ncfg_total,10);

[~,fov_coordinates_physical] = get_fov_coordinates(roixy,theta_local(1,:),theta_local(2,:),params);
fov_coordinates_physical_nm = 1e3*fov_coordinates_physical;
localizations_with_outliers(:,1) = ID;
localizations_with_outliers(:,2:3) = fov_coordinates_physical_nm(1:2,:)';
localizations_with_outliers(:,5) = framelist';
localizations_with_outliers(:,6) = CRLB_local(1,:)';
localizations_with_outliers(:,7) = CRLB_local(2,:)';
if contains(params.fitmodel,'xyz')
    localizations_with_outliers(:,4) = theta_local(3,:)'; 
    localizations_with_outliers(:,8) = CRLB_local(3,:)';
    localizations_with_outliers(:,9) = theta_local(4,:)';
    localizations_with_outliers(:,10) = theta_local(5,:)';
else
    localizations_with_outliers(:,4) = zeros(Ncfg_total,1);
    localizations_with_outliers(:,8) = zeros(Ncfg_total,1);
    localizations_with_outliers(:,9) = theta_local(3,:)';
    localizations_with_outliers(:,10) = theta_local(4,:)';
end

no_outliers = setdiff(1:Ncfg_total,outliers);
localizations = localizations_with_outliers(no_outliers,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.check_test
    disp(size(pos_store))
    figure;
    hold on;
    Legend = cell(1, Ncfg);
    for i = 1:Ncfg
        plot(1:max_local_iterations+1, pos_store(:, 1, i)');  % X
        Legend{i}=strcat('Spot #', num2str(i));
    end
    legend(Legend, 'Location', 'eastoutside');
    
    
    hold off;
    xlabel('Iteration');
    ylabel('X position (nm)');
    title('X Displacement over Iterations');

    figure;
    hold on;
    for i = 1:Ncfg
        plot(1:max_local_iterations+1, pos_store(:, 2, i)');  % Y
        Legend{i}=strcat('Spot #', num2str(i));
    end
    legend(Legend, 'Location', 'eastoutside');
    
    
    hold off;
    xlabel('Iteration');
    ylabel('Y position (nm)');
    title('Y Displacement over Iterations');

        figure;
    hold on;
    for ii = 1:Ncfg
        plot(1:max_local_iterations+1, pos_store(:, 3, ii)');  % Z
        Legend{i}=strcat('Spot #', num2str(i));
    end
    legend(Legend, 'Location', 'eastoutside');
    
    hold off;
    xlabel('Iteration');
    ylabel('Z position (nm)');
    title('Z Displacement over Iterations');

    %% observe individually
    % for ii = 1:Ncfg
    %     figure;
    %     imagesc(allspots(:,:,1,ii))
    %     add = sprintf('Spot # %i',ii);
    %     title(add, 'FontSize', 14);
    %     mu_spot = PSF_store(:,:,:,ii);
    %     dipshow(mu_spot)
    %     add_mu = sprintf('mu # %i',ii);
    %     title('mu # %i',ii)
    %     colormap parula
    %     diptruesize(2000)
    % end


    %% observe together
    row1 = [];
    row2 = [];
    row3 = [];
    row4 = [];
    t = tiledlayout(4,5);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    for ii = 1:Ncfg
        nexttile
        imagesc(allspots(:,:,1,ii))
      
       axis off
        mu_spot = PSF_store(:,:,:,ii);
        if ii <= 5

            row1 = cat(2,row1,mu_spot);

        end
        if ii > 5 && ii<=10

            row2 = cat(2,row2,mu_spot);

        end
        if ii > 10 && ii <= 15
            row3 = cat(2,row3,mu_spot);
          
        end
        if ii > 15 && ii <= 20
            row4 = cat(2,row4,mu_spot);
          
        end



    end
    grid = [row1; row2; row3; row4];
        size(grid)
        dipshow(grid(:,:,:),'lin','name','Measured ROIs vs. modeled PSFs (without outliers)');
        diptruesize(40000/params.Mx/2)
        colormap parula





    figure;
    hold on;
    for i = 1:Ncfg
        plot(1:max_local_iterations+1, grad_store(i, :, 3)');  % Gradient
        Legend{i}=strcat('Spot #', num2str(i));
    end
    legend(Legend, 'Location', 'eastoutside');
    hold off  
    xlabel('Iteration');
    ylabel('gradient Z');
    title('Z-Gradient per iteration');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % end of function

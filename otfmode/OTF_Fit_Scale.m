% Fit all ROIs with (or without) the estimated aberrations.
% The input data should be segmented
% You can choose between GPU, CPU or Matlab fitting


%%%%%%%%%%%%%%%%%%%%
% note: to run the test, you need to replace the parfor loop in
% 'local_update_otfmode' by a normal for loop (line 57), and comment the matrices
% storing the iterations in the while loop (lines 187-191)
%%%%%%%%%%%%%%%%%%%%%%



clear all
close all

input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\input data\transfer_3055388_files_4e4cf1e2\";
%input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\5_Localize";

path_aberrations = 'data/2_estimate_aberrations/estimated_aberrations_rep_0001_init_0001_segmentation_microtubules_3D_data0001.mat';

output_path = 'C:\Users\Naam\Documents\BEP\BEP Code Roman\data\Localize\';

params = set_parameters_microtubules_3D;
params.cpp = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.FlagOTF = true;

if params.FlagOTF
    params.fitmodel = 'xyz';
    path_aberrations = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_40.mat";%OTF

    params.use_fitted_aberrations = false; % If false, aberrations are set to 0, except for constant astigmatism for 3D data
    params.path_aberrations = path_aberrations;


    
    load(path_aberrations,'OTF_3d','PSFmodelmatrix','Xpatch','Ypatch','Zpatch','intensity');
    intensity = mean(intensity, [2, 3]);
    params.mult_run = true;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spots_per_batch = 2e4;

if params.cpp
    if strcmp(params.fitmodel,'xy')
        params.fitmodel = 'xy_constaberrations';
    elseif strcmp(params.fitmodel,'xyz')
        params.fitmodel = 'xyz_constaberrations';
    end
end

all_input_files = dir(input_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove current and parent directories
all_input_files = all_input_files(~[all_input_files.isdir]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr_of_files = size(all_input_files,1);
params.check_test = false; % use part of the dataset only to check your results.
% In order to be used, select one file, turn the parfor loop in
% local_update(_otfmode) into a for loop, and uncomment the arrays storing
% the gradients/spots


fprintf('Found %i files\n',nr_of_files);

tic_total = tic;

if params.check_test
    % picks one file 
    nr_of_files = randi([min(nr_of_files) max(nr_of_files)],1,1); % select a random file
end


for file_i= 78:100

    
    % Load data
    filename = all_input_files(file_i).name; 
    foldername = all_input_files(file_i).folder;
    path_input_data_full = fullfile(foldername, filename);    
    fprintf('Input segmented data = %s\n',path_input_data_full);

    loaded_data = load(path_input_data_full);
    allspots = loaded_data.allspots;
    roixy = loaded_data.roixy;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if params.mult_run && params.FlagOTF % this statement picks the OTF fit

        framelist = loaded_data.localizations(:,5);
        ID = loaded_data.localizations(:,1);
        remove = loaded_data.outliers;
        theta_fit = loaded_data.theta_update;
        
        locindex_file = 1:size(allspots,4);
        masklocs = ~ismember(locindex_file,remove);
        %% Makes Sure all the used variables are properly Seperated!
        %% Allspots represents all the spots of the original fit: you could sort it by the ID of the spot!
        Spots_fit = allspots(:,:,:,masklocs); 
        roixy_fit = roixy(:,masklocs); % We collect roixy_fit to describe all ROIs for the second fit, instead of for the first and second fit.
        theta_fit = theta_fit(:,masklocs);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % filtering out spots with negative Nph: not physically accurate!
        locs = loaded_data.localizations;
        no_negs = find(locs(:,9)>0);
        locs = locs(no_negs,:);
        Spots_fit = Spots_fit(:,:,:,no_negs);
        roixy_fit = roixy_fit(:,no_negs);
        theta_fit = theta_fit(:,no_negs);
        ID = ID(no_negs);
        framelist = framelist(no_negs);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % We are only interested in the spots within the z-limits of our
        % PSFs: filter out all non-compliant data
        %
        z_fit = [min(Zpatch) max(Zpatch)];
        sep = find((locs(:,4)>= z_fit(1)) & (locs(:,4)<= z_fit(2)));
        ID = ID(sep);
        framelist = framelist(sep);
        Spots_fit = Spots_fit(:,:,:,sep);
        roixy_fit = roixy_fit(:,sep);
        theta_fit = theta_fit(:,sep);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% This selects a random amount of spots from the file (20 spots from the first 1000
        if params.check_test
            % indx = randperm(1000,20);
            indx = 1:20;
            ID = ID(indx);
            framelist = framelist(indx);
            Spots_fit = Spots_fit(:,:,:,indx);
            roixy_fit = roixy_fit(:,indx);
            theta_fit = theta_fit(:,indx);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        params.select_spots = false;
        Ncfg_total = size(Spots_fit,4); %% Change to use only part of the data: for checking of results
        params.Ncfg_total = Ncfg_total;
        params.Ncfg_global = Ncfg_total;
        params.Ncfg = Ncfg_total;
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % We'll need to initialize some parameters, to properly use the
        % OTFmode

        if params.check_test
            spots_per_batch = 1e4;
            Ncfg_total = size(indx,2); %% Change to use only part of the data: for checking of results
            params.Ncfg_total = Ncfg_total;
            params.Ncfg_global = Ncfg_total;
            params.Ncfg = Ncfg_total;

        end

        % Much of this code is to ensure proper sizes
        Npatch = size(PSFmodelmatrix,5);%params.Npatch; just to ensure proper size
        params.Xpatch = Xpatch;
        params.Ypatch = Ypatch;
        params.Zpatch = Zpatch;
        params.Npatch =Npatch;
        % The initial parameters script needs to include some new
        % parameters:
        %%This can be removed
        xemit = 0;%params.xemit;
        yemit = 0;%params.yemit;
        zemit = 0;%params.zemit;
        params.x_emit = xemit;
        params.y_emit = yemit;
        params.z_emit = zemit;
        %% This is kept since the parameters file does define certain OTF size, but you might use others
        params.Notf = size(OTF_3d,1); %% Making sure parameters are appropriate size
        params.Notfz = size(OTF_3d,3);
        %%
        params.Mroix = params.Mx; %used in chirp-z transform
        params.Mroiy = params.My; %used in chirp-z transform
        params.roisamplingdistance = params.pixelsize; % sampling distance in ROI
        params.xroirange = params.roisamplingdistance*params.Mroix/2; % 1/2 of size PSF space in x (ROI)
        params.yroirange = params.roisamplingdistance*params.Mroiy/2; % 1/2 of size PSF space in y (ROI)
        params.fitmodel = 'xyz'; %ensures 3-D fitting
        %% The Parameter is originally zstage, which defaults to 0: Use the given emitter position from localization for an initial position estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        params.ztype = 'medium'; % ztype is initially stage, which sets all initial z-values to 0 for fitting
        % ztype = medium ensures the emission position is used as an
        % initial fit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        params.Mpsfx = size(PSFmodelmatrix,1);
        params.Mpsfy = size(PSFmodelmatrix,2);
        params.Mpsfz = size(PSFmodelmatrix,3);
        params.pos = linspace(-params.Mpsfz/2,params.Mpsfz/2,params.Mpsfz);
        params.zstage = zemit;
        params.compders =1;
        params.K=1; % z-stack = 1 (we only image projections of the PSF for now
        params.zspread = [min(Zpatch) max(Zpatch)]; % Using the OTF, we are limited to the spread of the 3-D PSF 
        % this limits all localization positions to be within these liits
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
    framelist = loaded_data.framelist;
    ID = loaded_data.ID;

    params.select_spots = false;
    Ncfg_total = size(allspots,4);
    params.Ncfg_total = Ncfg_total;
    params.Ncfg_global = Ncfg_total;
    params.Ncfg = Ncfg_total;
    end


    fprintf('\nPerform localization with %i ROIs\n',Ncfg_total);
    
    numparams = params.numparams;

    if params.use_fitted_aberrations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OTFmode doesn't need theta_glovbal
        fitted_aberrations_data = load(path_aberrations,'theta_full','max_merit_index');
        theta_global = fitted_aberrations_data.theta_full.global;

        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else    
        theta_global = zeros(params.numgammas,1);  
        theta_global(6) = 66; % Microtubules 3D average value
    end


    Nbatch = ceil(Ncfg_total/spots_per_batch);


    theta_init = zeros(params.numparams,Ncfg_total);
    outliers = [];
    mu = [];
    merit = zeros(Ncfg_total,1);
    meritoffset = zeros(Ncfg_total,1);
    Niters = [];


    % %% the previously stored values are re-used to again
    % %% localize! No need to use initialvalues_phasor!
    % % %preallocation of mu/Niters, and ensuring correct PSF size
    % % mu = zeros(params.Mx,params.My,params.K,Ncfg_total);
    % % Niters = zeros(Ncfg_total,1);
    params.Mpsfx = size(PSFmodelmatrix,1);
    params.Mpsfy = size(PSFmodelmatrix,2);
    params.Mpsfz = size(PSFmodelmatrix,3);



    % ID x y z framenumber CRLBx CRLBy CRLBz Nph bg
    localizations_with_outliers = zeros(Ncfg_total,10);
    localizations_with_outliers(:,1) = ID;
    localizations_with_outliers(:,5) = framelist;

    % Precalculation of OTFs
    if params.cpp && params.FlagOTF
        fitter_all = vectorfitter(allspots,roixy,params,params.cpp_fitmode);
        [otfs,spot_otf_indices] = fitter_all.update_otfs(theta_global,params.aberrations);
        delete(fitter_all)
    end
    
    % Fit in batches
    tic_file = tic;
    for batch = 1:Nbatch
        batch
        start = (batch-1)*spots_per_batch;
        stop = min(start+spots_per_batch,Ncfg_total);
        Nspots = stop - start;

        spots_batch = Spots_fit(:,:,:,start+1:stop);
        roixy_batch = roixy_fit(:,start+1:stop);
        

        if params.cpp
            fitter = vectorfitter(spots_batch,roixy_batch, params,params.cpp_fitmode);
            if params.FlagOTF
                spot_otf_indices_batch = spot_otf_indices(start+1:stop);
                fitter.set_otfs(otfs,spot_otf_indices_batch);
            else
                cpp_zernike_coeffs = fitter.update_zernike_coefficients(theta_global, params.aberrations);
            end
            theta_init_batch = fitter.theta_init();
            [theta_update_batch, Niters_batch] = fitter.local_update();
            dmudtheta_batch = fitter.dmudtheta;
            outliers_batch = fitter.outliers();
            crlb = fitter.fisher_crlb();
            mu_batch = fitter.mu();
            mu = cat(4,mu,mu_batch); 
            delete(fitter)
        else
            Ncfg = size(roixy_batch,2)
            params.Ncfg = Ncfg;
            params.Ncfg_total = Ncfg;
            thetastore_local = zeros(numparams,Ncfg,0);        
            meritstore_local = zeros(Ncfg,0);   
            alambda_local = ones(Ncfg,1)*params.alambda0_local;         
            Niters_batch = zeros(Ncfg,0);   
            iiter_total = 1;
            params.perform_final_iterations = false;
            

            flip_z_net = false;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if params.FlagOTF


                theta_init_batch = theta_fit(:,start+1:stop);
                %% For the initial fit, you don't yet have an otf: Let's keep these two similar functions separate
                [theta_update_batch,thetastore_local,localizations_temp,localizations_with_outliers_temp,meritstore_local,alambda_local,mu_batch,dmudtheta_batch,Niters_batch,outliers_batch] = ...
                local_update_otfmode2(theta_init_batch,thetastore_local,meritstore_local,alambda_local,spots_batch,roixy_batch,iiter_total,Niters_batch,flip_z_net,framelist(start+1:stop),ID(start+1:stop),OTF_3d,intensity,params); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else


                theta_init_batch = initialvalues_phasor(spots_batch,roixy_batch,theta_global,params);
                [theta_update_batch,thetastore_local,localizations_temp,localizations_with_outliers_temp,meritstore_local,alambda_local,mu_batch,dmudtheta_batch,Niters_batch,outliers_batch] = ...
                local_update(theta_init_batch,thetastore_local,theta_global,meritstore_local,alambda_local,spots_batch,roixy_batch,iiter_total,Niters_batch,flip_z_net,framelist(start+1:stop),ID(start+1:stop),params); 
            end
            crlb = localizations_with_outliers_temp(:,6:8)';
            outliers_batch = outliers_batch';
                
        end

        % Calculate the merit
        Ncfg = size(theta_update_batch,2);
        for jcfg = 1:Ncfg
            merit(start+jcfg) = likelihood(params,allspots(:,:,:,start+jcfg),mu_batch(:,:,:,jcfg),dmudtheta_batch(:,:,:,:,jcfg),params.varfit);
        end
        meritoffset(start+1:stop) = meritoffsetcalc(allspots(:,:,:,start+1:stop),params.varfit);

        no_outliers_temp = setdiff(1:Nspots,outliers_batch);
        outliers = cat(1,outliers,outliers_batch+start);
        theta_init(:,start+1:stop) = theta_init_batch;

        [~,xy_um] = get_fov_coordinates(roixy_batch,theta_update_batch(1,:),theta_update_batch(2,:),params);
        
        localizations_with_outliers(start+1:stop,2:3) = 1e3*xy_um';
        if contains(params.fitmodel,'xyz')
            localizations_with_outliers(start+1:stop,4) = theta_update_batch(3,:);
            localizations_with_outliers(start+1:stop,9:10) = theta_update_batch(4:5,:)';
            localizations_with_outliers(start+1:stop,6:8) = crlb(1:3,:)';
        else %xy
            localizations_with_outliers(start+1:stop,4) = zeros(Nspots,1);
            localizations_with_outliers(start+1:stop,9:10) = theta_update_batch(3:4,:)';
            localizations_with_outliers(start+1:stop,6:7) = crlb(1:2,:)';
            localizations_with_outliers(start+1:stop,8) = zeros(Nspots,1);
        end
        
        theta_update(:,start+1:stop) = theta_update_batch;
        theta_init(:,start+1:stop) = theta_init_batch;
        mu(:,:,:,start+1:stop) = mu_batch;
        Niters(start+1:stop) = Niters_batch;
        roixy_filter(:,start+1:stop) = roixy_batch; 

    end

    merit_corrected = merit+meritoffset;
    no_outliers = setdiff(1:Ncfg_total,outliers);
    localizations = localizations_with_outliers(no_outliers,:);

    t_file = toc(tic_file)
    locs_per_second = Ncfg_total/t_file



    filename = all_input_files(file_i).name; 
    path_output_data_full = strcat(output_path,'localization_OTF_',filename);
    save(path_output_data_full,'theta_update','theta_global','localizations','localizations_with_outliers','theta_init','outliers','roixy_fit','merit','merit_corrected','params','Spots_fit','mu','-v7.3')
    %save(path_output_data_full,'theta_update','theta_global','localizations','localizations_with_outliers','theta_init','outliers','roixy','merit','merit_corrected','params','-v7.3')
    
    fprintf('Localization with %i ROIs finished\n',Ncfg_total)

end

t_total = toc(tic_total)

%% Render final image
if params.show_plots
    x1 = localizations(:,2);
    y1 = localizations(:,3);
    z1  = localizations(:,4);
    
    locations = [x1 y1 z1];
   
    ncolors = 19;
    colorrange = [-500,400];
    
    xlim = [0.0 1.0];
    ylim = [0.0 1.0];
    pixelsize = 50

    sigma = 2.0;
    clim = 0.95; % the brightest (1-clim) fraction of pixels gets clipped 
    figure;
    render_3D_func(locations,ncolors,colorrange,xlim,ylim,pixelsize,sigma,clim,params);
end


%% Combine localizations from multiple files into one file
% 
% input_folder_combine = strcat(output_path,'localization*');
% output_path_final = strcat(output_path,'localization_combined.mat');
% 
% all_input_files = dir(input_folder_combine);
% nr_of_files = size(all_input_files,1);
% fprintf('Found %i files\n',nr_of_files);
% 
% localizations = [];
% %localizations_with_outliers = [];
% %theta_init = [];
% %theta = [];
% %outliers = [];
% roixy = [];
% %allspots = [];
% %mu = [];
% Ncfg_total = 0;
% 
% for file_i=1:nr_of_files
% %for file_i=1:7
% 
%     % Load data
%     filename = all_input_files(file_i).name; 
%     foldername = all_input_files(file_i).folder;
%     path_input_data_full = fullfile(foldername, filename);    
%     fprintf('Input data = %s\n',path_input_data_full);
%     loaded_data = load(path_input_data_full);
% 
%     localizations_temp = loaded_data.localizations;
%     %localizations_with_outliers_temp = loaded_data.localizations_with_outliers;
%     %theta_init_temp = loaded_data.theta_init;
%     %outliers_temp = loaded_data.outliers;
%     roixy_temp = loaded_data.roixy;
%     %allspots_temp = loaded_data.allspots;
%     %mu_temp = loaded_data.mu;
%     %theta_temp = loaded_data.theta_update;
%     params = loaded_data.params;
% 
%     localizations = cat(1,localizations,localizations_temp);
%     %localizations_with_outliers = cat(1,localizations_with_outliers,localizations_with_outliers_temp);
%     %theta_init = cat(2,theta_init,theta_init_temp);
%     %outliers_temp_adjusted = outliers_temp + Ncfg_total;
%     %outliers = cat(1,outliers,outliers_temp_adjusted);
%     %Ncfg_total = Ncfg_total + size(roixy_temp,2);
%     roixy = cat(2,roixy,roixy_temp);
%     %allspots = cat(4,allspots,allspots_temp);
%     %mu = cat(4,allspots,allspots_temp);
%     %theta = cat(2,theta,theta_temp);
% end
% 
% save(output_path_final,'localizations','roixy','params')


close all;
clear all;

path_input_data = "C:\Users\Naam\Documents\BEP\Code\vectorfit-master\vectorfit-master\BEP Code Roman\Extra input data\transfer_3055388_files_4e4cf1e2\";

path_output_data = "C:\Users\Naam\Documents\BEP\Code\vectorfit-master\vectorfit-master\data\2_estimate_aberrations\";

params = set_parameters_microtubules_3D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.cpp = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each init is done with the same subset, but with a different random 
% initialization.
% The highest likelihood solution is chosen, to avoid a local maximum.
% Recommended values:
% For 2D data: 10
% For 3D data: 1
nr_of_inits = 1;

% Repeat estimation of aberrations
% Each rep is done with a different subset of ROIs.
nr_of_reps = 1;

% Set up parallel pool
if isempty(gcp('nocreate'))
    number_of_parallel_workers = params.nr_of_parallel_workers;
    parpool('threads', number_of_parallel_workers);
end

for rep = 1:nr_of_reps
    for init=1:nr_of_inits
    
        % Select random input file from multiple segmentation files
        all_files = dir(strcat(path_input_data,'*'));
        n_files = numel(all_files);
        i_file = randi(n_files);
        filename = all_files(i_file).name; 
        params.filename_input = filename;
        foldername = all_files(i_file).folder;
        path_input_data_full = fullfile(foldername, filename);
        
        %path_input_data_full = path_input_data;
        fprintf('Input segmented data = %s\n',path_input_data_full);
        
        % Fit data
        tic;
        fprintf('\nStart fitting field dependent aberrations\n');
        Fit_NAT_data_otfmode
        toc;
        
        [~,filename_save,~]=fileparts(path_input_data_full); 
        path_output_data_full = strcat(path_output_data,'estimated_aberrations__otfmode','rep_',sprintf('%04d',rep),'_init_',sprintf('%04d',init),'_',filename_save,'.mat');
        save(path_output_data_full,'theta_full','thetastore_full','localizations_full','localizations_with_outliers_full','mu_full','dmudtheta_full','allspots','roixy','thetainit','params','max_merit_index','outliers_full','ID','meritstore_full','init');
        
        if params.show_plots
            show_fitted_data
        end
    end
end

% delete parpool
delete(gcp('nocreate'));

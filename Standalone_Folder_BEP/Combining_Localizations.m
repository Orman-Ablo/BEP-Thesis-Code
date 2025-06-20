%combining the localizations of all files

clear all
close all

%Define input and output file locations
input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Extra input data\transfer_3055388_files_4e4cf1e2\";
output_folder = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data";

%start the counter for total time
Tstart = tic;

all_input_files = dir(input_path);
all_input_files = all_input_files(~[all_input_files.isdir]);
nr_of_datasets = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_datasets);

%% set up parallel pool
nr_of_parallel_workers = 8;
if isempty(gcp('nocreate'))
    parpool('Threads', nr_of_parallel_workers);
end

loc_old = [];

%% loop over all files
for dataset_i=1:nr_of_datasets  
    dataset_name = all_input_files(dataset_i).name; 
    foldername = all_input_files(dataset_i).folder;
    input_filename = fullfile(foldername, dataset_name);
    fprintf('Input file = %s\n',dataset_name);
    loaded_input = load(input_filename);
    params = loaded_input.params;
    localizations = loaded_input.localizations;
    loc_new = localizations;
    loc_temp = [loc_old; loc_new];
    loc_old = loc_temp;


end
localizations = loc_temp;
size(localizations)
minimum = min(localizations(:,2:4));
minimum
maximum = max(localizations(:,2:4));
maximum
clear loc_old loc_new
path_output_data_full = strcat(output_folder,'localization_segmentation_','roi17_thr15p5_combined');
save(path_output_data_full,'localizations','params','-v7.3')
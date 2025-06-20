folder = "C:\Users\Naam\Documents\BEP\data\second_fit\";
folder_original = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Extra input data\transfer_3055388_files_4e4cf1e2\";
output = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Combining Results\";

all_input_files = dir(folder);
% Remove current and parent directories (show up as files across the
% directory
all_input_files = all_input_files(~[all_input_files.isdir]);

nr_of_files = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_files);

localizations = [];
roixy = [];
allspots = [];
mu = [];

for dataset_i = 1:100
    dataset_name = all_input_files(dataset_i).name; 
    foldername = all_input_files(dataset_i).folder;
    input_filename = fullfile(foldername, dataset_name);
    fprintf('Input file = %s\n',dataset_name);
    loaded_input = load(input_filename);
    localizations_temp = loaded_input.localizations;
    mu_temp = loaded_input.mu;
    allspots_temp = loaded_input.allspots;
    roixy_temp = loaded_input.roixy_fit;
    outliers = loaded_input.outliers;
    locindex_file = 1:size(mu_temp,4);
    masklocs = ~ismember(locindex_file,outliers);
    mu_temp  = mu_temp(:,:,:,masklocs);
    allspots_temp =  allspots_temp(:,:,:,localizations_temp(:,1));
    roixy_temp = roixy_temp(:,masklocs);
    params = loaded_input.params;
    roixy = cat(2,roixy,roixy_temp);
    mu = cat(4,mu,mu_temp);
    allspots = cat(4,allspots,allspots_temp);
    localizations = cat(1,localizations,localizations_temp);
end


filename = strcat(output,'Outliers_filtered');
save(filename,'localizations','roixy','mu','params','allspots')

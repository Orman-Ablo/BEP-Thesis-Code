folder = "C:\Users\Naam\Documents\BEP\data\second_fit\";
%folder_original = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Extra input data\transfer_3055388_files_4e4cf1e2\";
output = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Combining Results\";



all_input_files = dir(folder);
% Remove current and parent directories (show up as files across the
% directory
all_input_files = all_input_files(~[all_input_files.isdir]);

nr_of_files = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_files);

ID = [];
allspots = [];
for dataset_i = 1:100
    dataset_name = all_input_files(dataset_i).name; 
    foldername = all_input_files(dataset_i).folder;
    input_filename = fullfile(foldername, dataset_name);
    fprintf('Input file = %s\n',dataset_name);
    loaded_input = load(input_filename);
    localizations = loaded_input.localizations;
    params = loaded_input.params;
    allspots_temp = loaded_input.allspots;
    allspots = cat(4,allspots,allspots_temp);
    ID_temp = localizations(:,1);
    ID = [ID; ID_temp];
end

allspots = allspots(:,ID);
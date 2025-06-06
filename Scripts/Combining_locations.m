%folder_original = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Extra input data\transfer_3055388_files_4e4cf1e2\";
output = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Combining Results\";
folder = "C:\Users\Naam\Documents\BEP\data\final\";

all_input_files = dir(folder);
% Remove current and parent directories (show up as files across the
% directory
all_input_files = all_input_files(~[all_input_files.isdir]);

nr_of_files = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_files);

xstore = [];
ystore =[];
zstore = [];
for dataset_i = 1:nr_of_files
    dataset_name = all_input_files(dataset_i).name; 
    foldername = all_input_files(dataset_i).folder;
    input_filename = fullfile(foldername, dataset_name);
    fprintf('Input file = %s\n',dataset_name);
    loaded_input = load(input_filename);
    localizations = loaded_input.localizations;
    params = loaded_input.params;
    xtemp = localizations(:,2);
    ytemp = localizations(:,3);
    ztemp  = localizations(:,4);
    xstore = [xstore; xtemp];
    ystore = [ystore; ytemp];
    zstore = [zstore; ztemp];
end
locations = [xstore ystore zstore];
 
filename = strcat(output,'combined_locations_second_fit_final');
save(filename,'locations','params')
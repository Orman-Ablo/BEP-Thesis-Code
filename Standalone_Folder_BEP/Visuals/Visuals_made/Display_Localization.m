
%folder = "C:\Users\Naam\Documents\BEP\data\second_fit\";
folder = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Extra input data\transfer_3055388_files_4e4cf1e2";
%folder = "C:\Users\Naam\Documents\BEP\data\final\";

all_input_files = dir(folder);

% Remove current and parent directories (show up as files across the
% directory
all_input_files = all_input_files(~[all_input_files.isdir]);

nr_of_files = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_files);


xstore = [];
ystore =[];
zstore = [];
for dataset_i = 1
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


% locations = [x y z];
locations = [xstore ystore zstore];

ncolors = 19;
colorrange = [-500,400];

xlim = [0.0 1.0];
ylim = [0.0 1.0];
pixelsize = 50

sigma = 2.0;
clim = 0.95; % the brightest (1-clim) fraction of pixels gets clipped 
figure;
render_3D_func(locations,ncolors,colorrange,xlim,ylim,pixelsize,sigma,clim,params);
colorbar('southoutside')



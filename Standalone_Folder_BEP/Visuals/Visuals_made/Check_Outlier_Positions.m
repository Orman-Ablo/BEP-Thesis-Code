%determine the initial positions of the Outlier localizations


input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\data\Localize\";
input_path = "C:\Users\Naam\Documents\BEP\BEP Code Roman\data\Final_Localization";
folder_original = "C:\Users\Naam\Documents\BEP\BEP Code Roman\input data\transfer_3055388_files_4e4cf1e2\";
Npatch = 4;
check_file = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_40.mat";
load(check_file,'Xpatch','Ypatch','Zpatch')

all_input_files = dir(folder_original);
% Remove current and parent directories (show up as files across the
% directory
all_input_files = all_input_files(~[all_input_files.isdir]);

all_input_files_2 = dir(input_path);
% Remove current and parent directories (show up as files across the
% directory
all_input_files_2 = all_input_files_2(~[all_input_files_2.isdir]);

nr_of_files = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_files);
% xstore1 = [];
% ystore1 =[];
% zstore1 = [];
% roixy1 = [];
outlier = [];
firstfit = [];
secondfit  = [];

ratio1 = 0;
ratio2 = 0;
for dataset_i = 1:100
    dataset_name = all_input_files(dataset_i).name; 
    foldername = all_input_files(dataset_i).folder;
    input_filename = fullfile(foldername, dataset_name);
    fprintf('Input file = %s\n',dataset_name);
    loaded_input = load(input_filename);
    % localizations = loaded_input.localizations;
    localizations_with_outliers = loaded_input.localizations_with_outliers;
    outliers = loaded_input.outliers;

    ratio1 = ratio1 + (size(outliers)/size(localizations_with_outliers,1));

    first_fit = localizations_with_outliers(outliers,2:4);

    %% the second fit of data
    dataset_name = all_input_files_2(dataset_i).name; 
    foldername = all_input_files_2(dataset_i).folder;
    input_filename = fullfile(foldername, dataset_name);
    fprintf('Input file = %s\n',dataset_name);
    loaded_input = load(input_filename);
    %% Let's compare the outlier positions at the end of the first fit and second fit
    localizations_with_outliers_2 = loaded_input.localizations_with_outliers;
    outliers = loaded_input.outliers;

    ratio2 = ratio2 + (size(outliers)/size(localizations_with_outliers_2,1));
    second_fit  = localizations_with_outliers_2(outliers,2:4);


    % [match, fit] = ismember(localizations_with_outliers_2(:,1),outlier_locations(:,1));
    % first_fit_index =find(fit);
    % first_fit = localizations_with_outliers_2(first_fit_index,1:4);

    
    firstfit = [firstfit; first_fit];
    secondfit = [secondfit; second_fit];
end

disp(ratio1/100)
disp(ratio2/100)
out_fig = tiledlayout(Npatch,Npatch, 'TileSpacing','compact','Padding','tight');
axes_arr = gobjects(Npatch,Npatch);
legendHandles = [];
for ix = 1:Npatch
    for iy = 1:Npatch
         positionmask1 = (firstfit(:,1)>=Xpatch(ix))&(firstfit(:,1)<=Xpatch(ix+1))&...
                       (firstfit(:,2)>=Ypatch(iy))&(firstfit(:,2)<=Ypatch(iy+1));
         positionmask2 = (secondfit(:,1)>=Xpatch(ix))&(secondfit(:,1)<=Xpatch(ix+1))&...
                       (secondfit(:,2)>=Ypatch(iy))&(secondfit(:,2)<=Ypatch(iy+1));


         outlier_z1 = firstfit(positionmask1,3);
         outlier_z2 = secondfit(positionmask2,3);
         % max_z = max(Zpatch);
         % mean_z = mean(outlier_z);
         % min_z = min(Zpatch);
         % z_store(ix,iy,1) = min_z;
         % z_store(ix,iy,2) = max_z;
         % sprintf('for patch (%i,%i), the minimum is %i, and the maximum is %i',ix,iy,min_z,max_z)
         ax = nexttile;
         hold on
         hist1 = histogram(outlier_z1,53,'FaceAlpha', 0.3,'EdgeColor','none');
         hist2 = histogram(outlier_z2,53,'FaceAlpha',0.7,'EdgeColor','none');
         axes_arr(ix,iy) = ax;
         title(sprintf('Patch (%i,%i)',ix,iy))
         % xlabel('Z-pos(nm)')
         hold off
         if iy ~=1
             ax.YTickLabel = [];
         end

         if ix ~=4
             ax.XTickLabel = [];
         end

         
    end
end
for ii = 1:Npatch
    linkaxes(axes_arr(ii,:),'y')
end
xlabel(out_fig,'z-Positions (nm)')
ylabel(out_fig,'Counts')
sgtitle('Outlier Positions')
Lgnd = legend(legendHandles, 'Initial fit', 'OTF fit', 'Location', 'bestoutside');
fontsize(14,'points')
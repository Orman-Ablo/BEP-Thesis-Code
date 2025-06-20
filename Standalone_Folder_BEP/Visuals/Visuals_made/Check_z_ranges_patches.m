%% Now the files look very similar: points to a good fitting routine (or at least not bad)


input2 = "C:\Users\Naam\Documents\BEP\BEP Code Roman\data\Combined\localization_combined_second_fit_scaled.mat";
input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\data\Combined\localization_combined_initial_fit.mat";
load(input,'localizations_drift_corrected')
locations1 = localizations_drift_corrected(:,2:4);

load(input2,'localizations_drift_corrected')
locations2 = localizations_drift_corrected(:,2:4);
z_store2 = zeros(Npatch,Npatch,2);

Npatch = 4;
z_store = zeros(Npatch,Npatch,2);
check_file = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_40.mat";
load(check_file,'Xpatch','Ypatch','Zpatch')
zz = linspace(-1000,1000,50); % adjust to sample smoother/rougher
legendHandles = [];
figure;
Hist_fig = tiledlayout(Npatch,Npatch, 'TileSpacing','compact','Padding','compact');
axes_arr = gobjects(Npatch,Npatch);
for ix = 1:Npatch
    for iy = 1:Npatch
         positionmask1 = (locations1(:,1)>=Xpatch(ix))&(locations1(:,1)<=Xpatch(ix+1))&...
                       (locations1(:,2)>=Ypatch(iy))&(locations1(:,2)<=Ypatch(iy+1));
         locations_z1 = locations1(positionmask1,3);
         positionmask2 = (locations2(:,1)>=Xpatch(ix))&(locations2(:,1)<=Xpatch(ix+1))&...
                       (locations2(:,2)>=Ypatch(iy))&(locations2(:,2)<=Ypatch(iy+1));
         locations_z2 = locations2(positionmask2,3);
         % max_z = max(locations_z);
         % mean_z = mean(locations_z);
         % min_z = min(locations_z);
         % z_store(ix,iy,1) = min_z;
         % z_store(ix,iy,2) = max_z;
         % sprintf('for patch (%i,%i), the minimum is %i, and the maximum is %i',ix,iy,min_z,max_z)
         ax = nexttile;
         hold on
         histogram(locations_z1,zz,'FaceColor',[1 0 0],'FaceAlpha',0.5,'EdgeColor','none');
         histogram(locations_z2,zz,'FaceColor',[0 1 0],'FaceAlpha',0.5,'EdgeColor','none');
         hold off
         axes_arr(ix,iy) = ax;

         if iy ~=1
             ax.YTickLabel = [];
         end

         if ix ~=4
             ax.XTickLabel = [];
         end
         

         title(sprintf('Patch (%i,%i)',ix,iy))
         

         % xlabel('Z-pos(nm)')
        
         % legend('First fit','Second fit','Location','northwest')
         

    end
end
for ii = 1:Npatch
    linkaxes(axes_arr(ii,:),'y')
end
xlabel(Hist_fig,'z-Position (nm)')
ylabel(Hist_fig,'Counts')
Lgnd = legend(legendHandles, 'Initial fit', 'OTF fit', 'Location', 'bestoutside');
sgtitle('Histogram of z-Positions')
fontsize(14,'points')
% let's add in the same number of bins as in the PSFdata script

% figure;
% for ix = 1:Npatch
%     for iy = 1:Npatch
%          positionmask = (locations(:,1)>=Xpatch(ix))&(locations(:,1)<=Xpatch(ix+1))&...
%                        (locations(:,2)>=Ypatch(iy))&(locations(:,2)<=Ypatch(iy+1));
%          locations_z = locations(positionmask,3);
%          max_z = max(locations_z);
%          mean_z = mean(locations_z);
%          min_z = min(locations_z);
%          z_store2(ix,iy,1) = min_z;
%          z_store2(ix,iy,2) = max_z;
%          sprintf('for patch (%i,%i), the minimum is %i, and the maximum is %i',ix,iy,min_z,max_z)
%          subplot(Npatch,Npatch,(ix)+(iy-1)*4);
%          histogram(locations_z,zz)
%          title(sprintf('Patch (%i,%i)',ix,iy))
%          xlabel('Z-pos(nm)')
%          ylabel('spots')
% 
%     end
% end
% sgtitle('Z positions First Fit')
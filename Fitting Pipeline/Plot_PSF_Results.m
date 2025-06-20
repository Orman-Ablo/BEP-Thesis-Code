clear all 
close all
%input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\PSFdata_matrix_Shifted_Lidke3D.mat";
%Normalized PSF (after center shift)
input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\PSFdatamatrix_no_negative_Nph_Lidke3D.mat";
% input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\PSFdatamatrix_no_negative_Nph_x6_Lidke3D.mat";
% Normalized model PSF (after center shift)
%input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\PSFmodelmatrix_to_PSF_Lidke3D.mat";
%%
filename ="C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\localization_segmentation_roi17_thr15p5_combined.mat";

output_folder ='C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot through focus PSF, data and model side by side


load(filename,'localizations','params')
load(input,"PSFmodelmatrix","PSFdatamatrix","Xpatch","Ypatch","Zpatch",'chisquarematrix','numlocmatrix','PSF_Norm','bg_all','Nph_all')

show_psf = false;
show_otf = false;

%% keep in mind that now the PSFS and OTFs stored are a zstack: try to make it fit into the currently present code.

% parameters
pixelsize = params.pixelsize; % pixelsize in nm
upsfac = size(PSFdatamatrix,1)/params.Mx; % upsampling in xy
datacase = 'Lidke3D';

Npatch = size(PSFdatamatrix,5);
max(PSFdatamatrix,[],'all')

% takes the cumulative spots, subtracts background, divides photon count,
% averages over the number of spots
% Access the average PSF by removing the average background and photon
% count(?) for the spot
for ix = 1:Npatch
    for iy = 1:Npatch
        for iz = 1:size(PSFdatamatrix,3) 
            bg_all(iz,ix,iy) = bg_all(iz,ix,iy)/numlocmatrix(iz,ix,iy)/upsfac^2;
            Nph_all(iz,ix,iy) = Nph_all(iz,ix,iy);%/numlocmatrix(iz,ix,iy);
            PSFdatamatrix(:,:,iz,ix,iy) = (PSFdatamatrix(:,:,iz,ix,iy)-bg_all(iz,ix,iy)*numlocmatrix(iz,ix,iy))/ Nph_all(iz,ix,iy);% best performer: taking the mean of all results
        end

    end
end
    

PSF_final = PSFdatamatrix; % the weird shape is NOT due to the sum_shift_ups_PSF: it's got to be from either Nph or bg.


%% Normalize PSF?
L1 = sum((PSF_final), 'all');
PSF_final = PSF_final /L1;

% for ix = 1:Npatch
%     for iy = 1:Npatch
%         for iz = 1:size(PSFdatamatrix,3) 
%             PSF_final(:,:,iz,ix,iy) = PSF_final(:,:,iz,ix,iy)/sum(PSF_final(:,:,iz,ix,iy),'all');
%         end
% 
%     end
% end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define FOV patches

zz = (Zpatch(1:end-1)+Zpatch(2:end))/2;

[Nxy,~,Nz,Npatch,~] = size(PSFdatamatrix);
[Nxy_pix,~,~,~,~] = size(chisquarematrix);
upsfac = Nxy/Nxy_pix;
zz = (Zpatch(1:end-1)+Zpatch(2:end))/2;

% scale bar parameters
scalebarlength = 1;
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

locnumann = cell(Npatch,Npatch); %initialize cell array for annotation

% % position panels (x=rows=top-down, y=columns=left-right)
% panelposx = [0.71 0.49 0.27 0.05];
% panelposy = [0.03 0.26 0.49 0.73];
% panelsize = 130;

if show_psf
    % create movie object
    moviedir = 'C:\Users\Naam\Documents\BEP\BEP Code Roman\Figure Folder\';
    writerObjPSF = VideoWriter(strcat(moviedir,'throughfocusPSF_averaged',datacase,'.avi'));
    writerObjPSF.FrameRate = 1;
    open(writerObjPSF);
    
    for iz = 1:Nz
      fignum = iz;
      figure(fignum)
      figpos = [80 40 1440 720];
      set(gcf,'Position',figpos);
      tiledlayout(Npatch,Npatch,'TileSpacing','compact','Padding','compact')
      for ix = 1:Npatch
        for iy = 1:Npatch
          nexttile
          imagesc(PSF_final(:,:,iz,ix,iy))
          axis square
          axis off
          axis tight
          colormap parula
        
        end
      end
      sgtitle(strcat('z = ',num2str(zz(iz)),' nm'),'FontSize',20)
    %   if iz==1
        annotation('rectangle',[0.10 0.75 (0.096/2.13)*scalebarlength 0.01],'FaceColor','white','Color','white');
        annotation('textbox',[0.10 0.76 0.12 0.04],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
        for iy = 1:Npatch
          annotation('textbox',[0.08+(iy-1)*0.23 0.025 0.10 0.03],'String','data','FontSize',16,'Edgecolor','none','Color','black');
          annotation('textbox',[0.19+(iy-1)*0.23 0.025 0.10 0.03],'String','model','FontSize',16,'Edgecolor','none','Color','black');
        end
        for ix = 1:Npatch
          for iy = 1:Npatch
            annotation('textbox',[0.145+(iy-1)*0.229 0.945-0.23*(ix-1) 0.07 0.03],'String',strcat('(',num2str(ix),',',num2str(iy),')'),'FontSize',12,'Edgecolor','none','Color','black');
          end
        end
    %   end
    
      frame = getframe(gcf);
      writeVideo(writerObjPSF,frame);
    end
    figure;
    movie(frame)
    close(writerObjPSF);
    close all
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Figures in Slices:
%% Setting Grid Parameters as described
% this ensures proper use of parameter sizes across the file!
params.Mpsfx = size(PSF_final,1);%mPSFx;
params.Mpsfy = size(PSF_final,2);%mPSFy;
params.Mpsfz = size(PSF_final,3);% mPSFz;
mPSFx = params.Mpsfx;
mPSFy = params.Mpsfy;
mPSFz = params.Mpsfz;
nOTF = 48;
nOTFz = 40;


params.Notfx = nOTF;
params.Notfz = nOTFz;

x_rng = linspace(-(params.pixelsize*params.Mx)/2,(params.pixelsize*params.Mx)/2,mPSFx);
y_rng = linspace(-(params.pixelsize*params.My)/2,(params.pixelsize*params.My)/2,mPSFy);
z_rng = zz;% instead of giving intervals, this gives the midpoints of the intervals!
params.z_range = Zpatch(end);


% compders = 0,
%% Set of parameters not included in original 'Parameters' File, But needed for a 3D-OTF
%%%% this isn't accurate for the general files: it should match the size of
%%%% the ROIs!
params.psfsamplingdistance = 0.25*(params.lambda/params.NA/4); % sampling distance in PSF space
params.psfsamplingdistancez = 0.1*(params.lambda/params.NA/4); % sampling distance in PSF space
%%%%




params.roisamplingdistance = params.pixelsize; % sampling distance in ROI
params.psfxrange = params.psfsamplingdistance*mPSFx/2; % 1/2 of size PSF space in x
params.psfyrange = params.psfsamplingdistance*mPSFy/2; % 1/2 of size PSF space in y
params.psfzrange = params.psfsamplingdistancez*mPSFz/2;
params.Notf = nOTF; % #sampling points in OTF space, issue #1 form standard params
params.spat_freq_cutoff = params.NA*2/params.lambda; % Lateral cutoff frequency of parameters
n_imm = params.refimmnom;
params.spat_freq_cutoffz = (n_imm - sqrt(n_imm^2 - (params.NA)^2))/params.lambda;% axial cutoff frequency of parameters
params.spat_freq_res = 1/(params.pixelsize*params.Mx); % resolution based on the size of the ROIs
params.spat_freq_resz = 1/(Zpatch(2)-Zpatch(1));% z-resolution based on the ROIs
params.OTF_pixel_cutoff= ceil(params.spat_freq_cutoff/params.spat_freq_res); % smallest number of pixels that encapsulate the cutoff frequency
params.OTF_pixel_cutoffz= ceil(params.spat_freq_cutoffz/params.spat_freq_resz);% smallest number of pixels that encapsulate the cutoff frequency in the axial direction
params.OTF_cutoff = params.OTF_pixel_cutoff*nOTF/params.Mx; % scaled OTF cutoff based on final size
params.OTF_cutoffz = params.OTF_pixel_cutoffz*nOTFz/mPSFz; % scaled OTF cutoff based on final size

compders = 0;


%% collecting Arrays
OTF_3d = zeros(nOTF,nOTF,nOTFz,Npatch,Npatch);
OTF_3d_abs = zeros(nOTF,nOTF,nOTFz,Npatch,Npatch);

% % try to pad the PSF to preserve the shape
% NzPad = 10; % or even 20
% PSF_final = padarray(PSF_final, [0, 0, NzPad,0,0], 0, 'both');
% params.Mpsfz = (params.Mpsfz)+NzPad*2;

params.fitmodel = 'xyz';
for xi = 1:Npatch
    for yi = 1:Npatch
          OTF_3d(:,:,:,xi,yi) = get_otf(PSF_final(:,:,:,xi,yi),params);%calculate the 3D OTf patch-by-patch
          OTF_3d_abs(:,:,:,xi,yi)= abs(OTF_3d(:,:,:,xi,yi));
    end
end
% % restore the original PSF shape
% PSF_final  = PSF_final(:,:,NzPad+1:end-NzPad,:,:);
% size(PSF_final)


% Check energy conservation!
% cztfuncn2 maintains the total energy of the system!
% OTF = OTF_3d(:,:,:,2,2);
% E_OTF = sum(abs(OTF).^2,'all');
% disp(E_OTF)
% PSF = PSF_final(:,:,:,2,2);
% E_PSF = sum(abs(PSF).^2,'all');
% disp(E_PSF)
% disp(E_PSF/E_OTF)
% OTF = OTF*(sqrt(E_PSF/E_OTF));
% disp(sum(abs(OTF).^2,'all'))




if show_otf
    % Movie initialization 
    moviedir = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Figure Folder\";
    writerObjOTF = VideoWriter(strcat(moviedir,'throughfocusOTF_all_files_',num2str(nOTF),'_',num2str(nOTFz),'.avi'));
    writerObjOTF.FrameRate = 4;
    open(writerObjOTF);
    
    
    %% set limits of transfer functions 
    % clim_psf = [min(PSF_final, [], 'all') max(PSF_final,[],'all')];
    clim_psf = [0 max(PSF_final,[],'all')];
    clim_otf = [min(OTF_3d_abs,[],'all') max(OTF_3d_abs, [], 'all')];
    
    %movie file of the OTF
    qx_range_OTF=linspace(-params.spat_freq_cutoff,params.spat_freq_cutoff,nOTF+1);
    qy_range_OTF=linspace(-params.spat_freq_cutoff,params.spat_freq_cutoff,nOTF+1);
    z_range_otf =linspace(-params.spat_freq_cutoffz,params.spat_freq_cutoffz,nOTFz+1);
    for zi = 1:nOTFz
        figure
        figpos = [80 40 1080 720];
        set(gcf,'Position',figpos);
        for xi = 1:Npatch
            for yi = 1:Npatch
            subplot(Npatch,Npatch,Npatch*(xi-1)+yi)
            imagesc(qx_range_OTF,qy_range_OTF,OTF_3d_abs(:,:,zi,xi,yi),clim_otf)
            title(sprintf('Modeled |OTF|'))
            xlabel('qx')
            ylabel('qy')
            zlabel('qz')
            axis square
            %axis off
            ax = gca;
            ax.FontSize = 8;
    
            end
        end
        sgtitle(strcat('f =  ',num2str(z_range_otf(zi)),'  to  ',num2str(z_range_otf(zi+1)),'  cycles/nm'),'FontSize',20);
        frame(zi) = getframe(gcf);
        writeVideo(writerObjOTF,frame(zi));
        %figname = sprintf('OTF Full_Figure_%i.fig',zi);
        %savefig(gcf, fullfile("C:/Users/Naam/Documents/BEP/Code/vectorfit-master/vectorfit-master/BEP Code Roman/Figure Folder/", figname))
    end
    
    %Show the video
    figure;
    movie(frame)
    close(writerObjOTF);
    close all
end



%% Show the PSF and OTF in slices
% Helps us visualize the data in 3D, alongside the Throughfocus Movies
% size_spots = size(spots_count);


% [X_PSF, Y_PSF, Z_PSF] = meshgrid(1:size(mu_final, 1), 1:size(mu_final, 2), 1:size(mu_final, 3));
% [X_OTF, Y_OTF, Z_OTF] = meshgrid(1:size(OTF_3d_abs, 1), 1:size(OTF_3d_abs, 2), 1:size(OTF_3d_abs, 3));
qx_range_OTF=linspace(-params.spat_freq_cutoff,params.spat_freq_cutoff,nOTF);
qy_range_OTF=linspace(-params.spat_freq_cutoff,params.spat_freq_cutoff,nOTF);
qz_range_OTF =linspace(-params.spat_freq_cutoffz,params.spat_freq_cutoffz,nOTFz);


upsfac= params.Mpsfx/params.Mx;
%
 directory  ="C:\Users\Naam\Documents\BEP\BEP Code Roman\Figure Folder\";
 figure;
 for xi = 1:Npatch
     for yi = 1:Npatch
        subplot(Npatch,Npatch,Npatch*(xi-1)+yi)
        plot_matrix = PSF_final(:,:,:,xi,yi);
        %% Finds the index and slice postion that maps to the highest intensity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Max_PSF, Index_PSF] = max(plot_matrix,[],'all');
        [x_slc, y_slc, z_slc] = ind2sub(size(plot_matrix),find(plot_matrix == Max_PSF));
        xslice_PSF = ((x_slc*pixelsize)/upsfac)-(params.pixelsize*params.Mx)/2;
        yslice_PSF = ((y_slc*pixelsize/upsfac))-(params.pixelsize*params.Mx)/2;
        zslice_PSF = (z_slc*pixelsize/upsfac)-zz(end);      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plt = slice(x_rng,y_rng,z_rng,plot_matrix,xslice_PSF, yslice_PSF, zslice_PSF);
        set(plt,'edgecolor','none')
        grid off
        xlabel('x (nm)')
        ylabel('y (nm)')
        zlabel('z (nm)')
        title(sprintf('PSF: Patch (%d,%d)',xi,yi))
        
        %clim(clim_psf)
     end
 end
 file_PSF = fullfile(directory, sprintf('PSF_sliced_%d_%d.fig',mPSFx,mPSFz));
 savefig(file_PSF);
 figure;
 for xi = 1:Npatch
     for yi = 1:Npatch
        subplot(Npatch,Npatch,Npatch*(xi-1)+yi)
        plot_otf = OTF_3d_abs(:,:,:,xi,yi);
        [Max_OTF, Index_OTF] = max(plot_otf,[],'all');
        [qx_slc, qy_slc, qz_slc] = ind2sub(size(plot_otf),find(plot_otf == Max_OTF));
        qxslice_OTF = qx_slc*(qx_range_OTF(2)-qx_range_OTF(1))-qx_range_OTF(end);
        qyslice_OTF = qy_slc*(qy_range_OTF(2)-qy_range_OTF(1))-qy_range_OTF(end);
        qzslice_OTF = qz_slc*(qz_range_OTF(2)-qz_range_OTF(1))-qz_range_OTF(end);
        plt = slice(qx_range_OTF,qy_range_OTF,qz_range_OTF,plot_otf,qxslice_OTF, qyslice_OTF, qzslice_OTF);
        set(plt,'edgecolor','none')
        xlabel('qx')
        ylabel('qy')
        zlabel('qz')
        title(sprintf('OTF: Patch (%d,%d)',xi,yi))

        
        grid off
        %clim(clim_otf)
        %set(gca,'CLim',[minlev maxlev])
     end
 end
 file_OTF = fullfile(directory, sprintf('OTF_sliced_%d_%d.fig',nOTF,nOTFz));
 savefig(file_OTF);


 % Alternative: look at the xz and yz profiles seperately: easier to
 % visualize!
  figure;
  for xi = 1:Npatch
     for yi = 1:Npatch
         plot_matrix = PSF_final(:,:,:,xi,yi);
         subplot(Npatch,Npatch,Npatch*(xi-1)+yi)
         [Max_PSF, Index_PSF] = max(plot_matrix,[],'all');
         [x_slc, y_slc, z_slc] = ind2sub(size(plot_matrix),find(plot_matrix == Max_PSF));
         img = plot_matrix(:,y_slc,:);
         img=squeeze(img)';
         img = flip(img,1);
 
         zz = flip(zz);
         imagesc(y_rng,zz,img)
         xlabel('y (nm)')
         ylabel('z (nm)')
         title(sprintf('PSF (yz): Patch (%d,%d)',xi,yi))
     end
  end
  figure;
  for xi = 1:Npatch
     for yi = 1:Npatch
          plot_matrix = PSF_final(:,:,:,xi,yi);
         subplot(Npatch,Npatch,Npatch*(xi-1)+yi)
         [Max_PSF, Index_PSF] = max(plot_matrix,[],'all');
         [x_slc, y_slc, z_slc] = ind2sub(size(plot_matrix),find(plot_matrix == Max_PSF));
         img = plot_matrix(:,y_slc,:);
         img=squeeze(img)';
         img = flip(img,1);
     
         zz = flip(zz);
         imagesc(x_rng,zz,img)
         xlabel('x (nm)')
         ylabel('z (nm)')
         title(sprintf('PSF (xz): Patch (%d,%d)',xi,yi))
     end
  end

    figure;
  for xi = 1:Npatch
     for yi = 1:Npatch
         plot_otf = OTF_3d_abs(:,:,:,xi,yi);
         subplot(Npatch,Npatch,Npatch*(xi-1)+yi)
         [Max_OTF, Index_PSF] = max(plot_otf,[],'all');
         [x_slc, y_slc, z_slc] = ind2sub(size(plot_otf),find(plot_otf == Max_OTF,1));
         img = plot_otf(:,y_slc,:);
         img=squeeze(img)';
         img = flip(img,1);
         imagesc(qy_range_OTF,qz_range_OTF,img)
         xlabel('qy (nm)')
         ylabel('qz (nm)')
         title(sprintf('OTF (qy/qz): Patch (%d,%d)',xi,yi))
     end
  end



%% let's look at the intensity profile of the middle grid across z:
% it should be negatively parabolic over z: I(z) = (I_0*(1/(1+((z/zo)^2))
intensity = zeros(mPSFz,Npatch,Npatch);
figure;
for zi = 1:mPSFz
    for xi = 1:Npatch
        for yi = 1:Npatch

            intensity(zi,xi,yi) = sum(PSF_final(:,:,zi,xi,yi),'all');

        end
    end
end
X = linspace(0,mPSFz,mPSFz);
min_I =min(intensity,[],'all');
max_I =max(intensity,[],'all');
lims = [min_I max_I];
X = transpose(X);
for xi = 1:Npatch
    for yi = 1:Npatch
       intensity_iter = intensity(:,xi,yi);
       min_I =min(intensity_iter,[],'all');
       max_I =max(intensity_iter,[],'all');
       lims = [min_I max_I];
       subplot(Npatch,Npatch,Npatch*(xi-1)+yi)
       plot(zz,intensity_iter)
       ylim(lims)
       title(sprintf('Axial Intensity Profile (%d,%d)',xi,yi))
       xlabel('Axial Coordinate (nm)')
       ylabel('Intensity (I)')     
    end
end


path_output_data_full = strcat(output_folder,'OTF_PSF_Store_',num2str(nOTFz));
save(path_output_data_full,'PSF_final','PSFdatamatrix','PSFmodelmatrix','Xpatch','Ypatch','Zpatch','chisquarematrix','numlocmatrix','OTF_3d','OTF_3d_abs','params','Nph_all','bg_all','-v7.3')


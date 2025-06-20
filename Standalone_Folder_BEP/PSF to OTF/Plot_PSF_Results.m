clear all 
close all

input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\PSFdatamatrix_no_negative_Nph_Lidke3D.mat";
%input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\PSFdatamatrix_2muLidke3D.mat";
%input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\PSFdatamatrix_1200nm_Lidke3D.mat";
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
% count for each grid
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 = sum((PSF_final), 'all');
% PSF_final = PSF_final /L1;

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
pix_siz = Zpatch(2)-Zpatch(1);


% compders = 0,
%% Set of parameters not included in original 'Parameters' File, But needed for a 3D-OTF
params.roisamplingdistance = params.pixelsize; % sampling distance in ROI
params.psfxrange = pix_siz*mPSFx/2; % 1/2 of size PSF space in x
params.psfyrange = pix_siz*mPSFy/2; % 1/2 of size PSF space in y
params.psfzrange = pix_siz*mPSFz/2; % 1/2 of size PSF space in y
params.Notf = nOTF; % #sampling points in OTF space, issue #1 form standard params

n_imm = params.refimmnom;


params.spat_freq_res = 1/(params.pixelsize*params.Mx); % resolution based on the size of the ROIs
params.spat_freq_resz = 1/(pix_siz*mPSFz); %qz-resolution


axis_x_y = -floor(params.Notf/2): ceil(params.Notf/2)-1;
qx_qy = axis_x_y*params.spat_freq_res;

% [qx,qy]  = meshgrid(qx_qy,qx_qy);
qx = qx_qy;
qy = qx_qy;

axis_z = -floor(params.Notfz/2):ceil(params.Notfz/2)-1;
qz = axis_z*params.spat_freq_resz;


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
          
          
          % OTF = get_otf(PSF_final(:,:,:,xi,yi),params);%calculate the 3D OTf patch-by-patch
          % PSF = PSF_final(:,:,:,xi,yi);
          % energy_PSF = sum(abs(PSF(:)).^2);
          % energy_OTF = sum(abs(OTF(:)).^2);
          % OTF_3d(:,:,:,xi,yi) = OTF * sqrt(energy_PSF / energy_OTF);
          OTF_3d_abs(:,:,:,xi,yi)= abs(OTF_3d(:,:,:,xi,yi));
    end
end


% % restore the original PSF shape
% PSF_final  = PSF_final(:,:,NzPad+1:end-NzPad,:,:);
% size(PSF_final)

% PSF = PSF_final(:,:,:,1,1);
% OTF = OTF_3d(:,:,:,1,1);
% 
% energy_PSF = sum(abs(PSF(:)).^2);
% energy_OTF = sum(abs(OTF(:)).^2);
% disp(energy_PSF);
% disp(energy_OTF);


% OTF = OTF * sqrt(energy_PSF / energy_OTF);
% 
% energy_PSF_after = sum(abs(PSF(:)).^2);
% energy_OTF_after = sum(abs(OTF(:)).^2);
% disp(energy_PSF_after)
% disp(energy_OTF_after)


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

if show_psf
    % create movie object
    moviedir = 'C:\Users\Naam\Documents\BEP\BEP Code Roman\Visuals\Figure Folder\';
    writerObjPSF = VideoWriter(strcat(moviedir,'throughfocusPSF_averaged_2mu',datacase,'.avi'));
    writerObjPSF.FrameRate = 1;
    open(writerObjPSF);
    
    for iz = 1:Nz
      fignum = iz;
      figure(fignum)
      figpos = [80 40 1706.7 946];
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
        annotation('textbox',[0.10 0.753 0.12 0.04],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
        for iy = 1:Npatch
          annotation('textbox',[0.14+(iy-1)*0.23 0.025 0.10 0.03],'String','data','FontSize',16,'Edgecolor','none','Color','black');
        end
        for ix = 1:Npatch
          for iy = 1:Npatch
            annotation('textbox',[0.145+(iy-1)*0.23 0.955-0.235*(ix-1) 0.07 0.03],'String',strcat('(',num2str(ix),',',num2str(iy),')'),'FontSize',12,'Edgecolor','none','Color','black');
          end
        end

          for ix = 1:Npatch
            for iy = 1:Npatch
              numlocstring = strcat('K=',num2str(numlocmatrix(iz,ix,iy)));
              locnumann{ix,iy} = annotation('textbox',[0.10+(iy-1)*0.23 0.93-0.235*(ix-1) 0.07 0.03],'String',numlocstring,'FontSize',12,'Edgecolor','none','Color','white');
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



if show_otf
    % Movie initialization 
    moviedir = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Visuals\Figure Folder\";
    writerObjOTF = VideoWriter(strcat(moviedir,'throughfocusOTF_all_files_1200nm_',num2str(nOTF),'_',num2str(nOTFz),'.avi'));
    writerObjOTF.FrameRate = 4;
    open(writerObjOTF);
    
    
    %% set limits of transfer functions 
    % clim_psf = [min(PSF_final, [], 'all') max(PSF_final,[],'all')];
    clim_psf = [0 max(PSF_final,[],'all')];
    clim_otf = [min(OTF_3d_abs,[],'all') max(OTF_3d_abs, [], 'all')];
    
    
  
    scale = 0.1566;
    scalebarlength = 0.005; % in nm^-1
    scalebarstring = strcat(num2str(scalebarlength),' (1/nm)');
    pixel_bar = scalebarlength / params.spat_freq_res;
    bar_scale = (pixel_bar / nOTFz) * scale;

    
    %movie file of the OTF
    for zi = 1:nOTFz-1
        figure
        figpos = [80 40 1080 720];
        set(gcf,'Position',figpos);
        for xi = 1:Npatch
            for yi = 1:Npatch
            subplot(Npatch,Npatch,Npatch*(xi-1)+yi)
            imagesc(qx,qy,OTF_3d_abs(:,:,zi,xi,yi),clim_otf)
            title(sprintf('Modeled |OTF|'))
            xlabel('qx (nm^{-1})')
            ylabel('qy (nm^{-1})')
            zlabel('qz (nm^{-1})')
            axis square
            %axis off
            ax = gca;
            ax.FontSize = 8;
    
            end
        end
        sgtitle(strcat('f =  ',num2str(qz(zi)),'  to  ',num2str(qz(zi+1)),'  (nm^{-1})'),'FontSize',20);
        annotation('rectangle',[0.17 0.753 bar_scale 0.005],'FaceColor','white','Color','white');
        annotation('textbox',[0.165 0.7545 0.12 0.04],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
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
upsfac= params.Mpsfx/params.Mx;
%
 directory  ="C:\Users\Naam\Documents\BEP\BEP Code Roman\Visuals\Figure Folder\";
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
        qxslice_OTF = qx_slc*(qx(2)-qx(1))-qx(end);
        qyslice_OTF = qy_slc*(qy(2)-qy(1))-qy(end);
        qzslice_OTF = qz_slc*(qz(2)-qz(1))-qz(end);
        plt = slice(qx,qy,qz,plot_otf,qxslice_OTF, qyslice_OTF, qzslice_OTF);
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


scale = 0.1262;
scalebarlength = 500; % in nm
scalebarstring = strcat(num2str(scalebarlength),' (nm)');
pixel_bar = scalebarlength / (Zpatch(2)-Zpatch(1));
bar_scale = (pixel_bar / mPSFx) * scale;


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
         annotation('rectangle',[0.15 0.77 bar_scale 0.007],'FaceColor','white','Color','white');
         annotation('textbox',[0.15 0.767 0.12 0.04],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
         
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
         annotation('rectangle',[0.15 0.77 bar_scale 0.007],'FaceColor','white','Color','white');
         annotation('textbox',[0.15 0.767 0.12 0.04],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
        
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
         imagesc(qy,qz,img)
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
save(path_output_data_full,'PSF_final','PSFdatamatrix','PSFmodelmatrix','Xpatch','Ypatch','Zpatch','chisquarematrix','numlocmatrix','OTF_3d','OTF_3d_abs','params','Nph_all','bg_all','intensity','-v7.3')


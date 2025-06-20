
input = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\OTF_PSF_Store_40.mat";
load(input,'PSF_final','OTF_3d','params','Zpatch')

mPSFx = params.Mpsfx;
mPSFy = params.Mpsfy;
mPSFz = params.Mpsfz;
nOTF = params.Notfx;
nOTFz = params.Notfz;
[Nxy,~,Nx,Npatch,~] =size(PSF_final);
zz = (Zpatch(1:end-1)+Zpatch(2:end))/2;


show_psf = false;
show_otf = false;
slice_otf = false;
slice_psf = false;
project_psf = true;
project_otf = true;
intensity = true;

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

if slice_psf
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
        xslice_PSF = ((x_slc*params.pixelsize)/upsfac)-(params.pixelsize*params.Mx)/2;
        yslice_PSF = ((y_slc*params.pixelsize/upsfac))-(params.pixelsize*params.Mx)/2;
        zslice_PSF = (z_slc*params.pixelsize/upsfac)-zz(end);      
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

end

if slice_otf
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
end

if project_psf
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

         imagesc(x_rng,zz,img)
         xlabel('x (nm)')
         ylabel('z (nm)')
         title(sprintf('PSF (xz): Patch (%d,%d)',xi,yi))
         annotation('rectangle',[0.15 0.77 bar_scale 0.007],'FaceColor','white','Color','white');
         annotation('textbox',[0.15 0.767 0.12 0.04],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
        
     end
  end
end

if project_otf

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
end


if intensity
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
end
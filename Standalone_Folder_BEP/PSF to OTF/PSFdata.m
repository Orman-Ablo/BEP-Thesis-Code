% This script is for evaluating the 3D-PSF across the FOV from the
% localization data.
%
% Sjoerd Stallinga, TUD, 2019-2024

close all
clear all

%%
% read in localization data

fprintf('load localization data ...\n')

datacase = 'Lidke3D';
%datacase = 'Fu3D';

% read-in data and set parameters
switch datacase
  case 'Lidke3D'
    locdatadir = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\";
    % locdatafilename = 'roi17_thr15p5_combined.mat';
    locdatafilename = "roi17_thr15p5_combined.mat";
    numfiles = 100;
    Nxy = 17;
  case 'Fu3D'
    locdatadir = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu\data_analysis\analysis_vectorfit\localization\localization_with_aberrations\localization_combined\';
    locdatafilename = 'roi17_thr24_NUP96_SNP647_astigmatism_3D_1608_10ms_hama_mm_1800mW_3_MMStack_combined.mat';
    numfiles = 122;
    Nxy = 17;
end
filename = strcat(locdatadir,'localization_segmentation_',locdatafilename);

load(filename,'localizations','params')

% parameters
pixelsize = params.pixelsize; % pixelsize in nm
upsfac = 7; % upsampling in xy

% setup array for properly shifting the ROIs to the center
numlocs = size(localizations,1);
locindex = localizations(:,1);
coords = localizations(:,2:4);
clear localizations

% define FOV patches
zrange = 1000; % total axial range
Npatch = 4;
focalslice = pixelsize/upsfac;
Nz = round(zrange/focalslice);
minxyz = min(coords);
maxxyz = max(coords);
minxyz(3) = -Nz*pixelsize/upsfac/2;
maxxyz(3) = Nz*pixelsize/upsfac/2;

Xpatch = minxyz(1)+(0:Npatch)*(maxxyz(1)-minxyz(1))/Npatch;
Ypatch = minxyz(2)+(0:Npatch)*(maxxyz(2)-minxyz(2))/Npatch;
Zpatch = minxyz(3)+(0:Nz)*(maxxyz(3)-minxyz(3))/Nz;
zz = (Zpatch(1:end-1)+Zpatch(2:end))/2;

%%
% loop over ROI and PSF model files

% path to model PSF data
switch datacase
  case 'Lidke3D'
    modeldatadir = 'C:\Users\Naam\Documents\BEP\BEP Code Roman\input data\transfer_3055388_files_4e4cf1e2\';
    modeldatafilename = 'roi17_thr15p5_Data0'; 
  case 'Fu3D'
    % modeldatadir = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu\data_analysis\analysis_vectorfit\localization\localization_with_aberrations\localization_allfiles\';
    modeldatadir = 'U:\SMLMaberrations\NATaberrations\Code\PSF analysis\localization_allfiles\';
    modeldatafilename = 'roi17_thr24_NUP96_SNP647_astigmatism_3D_1608_10ms_hama_mm_1800mW_3_MMStack_Default_'; 
end

% selection of data files
allfiles = 1:numfiles; % maximum amount, this is all data
% allfiles = 12:12;

% define array for average PSF from data
PSFdatamatrix = zeros(upsfac*Nxy,upsfac*Nxy,Nz,Npatch,Npatch);
PSFmodelmatrix = zeros(upsfac*Nxy,upsfac*Nxy,Nz,Npatch,Npatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSF_Norm = zeros(upsfac*Nxy,upsfac*Nxy,Nz,Npatch,Npatch);

ratio_min_all = zeros(size(allfiles));
ratio_max_all = zeros(size(allfiles));
ratio_mean_all = zeros(size(allfiles));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numlocmatrix = zeros(Nz,Npatch,Npatch);
bg_all  = zeros(Nz,Npatch,Npatch);
Nph_all = zeros(Nz,Npatch,Npatch);

% define array for computing the chi-square goodness of fit measure
chisquarematrix = zeros(Nxy,Nxy,Nz,Npatch,Npatch);

for jfile = allfiles
  
  % load ROI and model data
  switch datacase
    case 'Lidke3D'
      kfile = jfile;
    case 'Fu3D'
      kfile = jfile-1;
  end
  if kfile>=100
    filenumberindicator = num2str(kfile);
  elseif kfile>=10
    filenumberindicator = strcat('0',num2str(kfile));
  else
    filenumberindicator = strcat('00',num2str(kfile));
  end
    
  % load PSF model data
  fprintf('load PSF data file %i...\n',jfile)
  switch datacase
    case 'Lidke3D'
      filename = strcat('localization_segmentation_',modeldatafilename,filenumberindicator,'.mat'); 
    case 'Fu3D'
      filename = strcat('localization_segmentation_',modeldatafilename,filenumberindicator,'.ome.mat'); 
  end
  filename = strcat(modeldatadir,filename);
  load(filename,'allspots','localizations','mu','outliers','roixy')
  allmodelPSFs = squeeze(mu);
  allspots = squeeze(allspots);
  roixy = roixy';


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % filtering out spots with negative Nph:
  no_negs = find(localizations(:,9)>0);
  bg = localizations(no_negs,10);
  Nph = localizations(no_negs,9);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  locindex_file = 1:size(allspots,3);
  masklocs = ~ismember(locindex_file,outliers);
  coords_file = localizations(:,2:4);
  allspots = allspots(:,:,masklocs);
  allmodelPSFs = allmodelPSFs(:,:,masklocs);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  allspots = allspots(:,:,no_negs);
  allmodelPSFs = allmodelPSFs(:,:,no_negs);
  coords_file = coords_file(no_negs,:);

  allspotsx = sum(allspots,1);
  allspotssum = sum(allspotsx,2);
  mux = sum(allmodelPSFs,1);
  musum = sum(mux,2);
  ratio = musum./allspotssum;
  % sprintf('max');
  % disp(max(ratio));% compare mu to allspots: see if you can use it!
  % sprintf('min');
  % disp(min(ratio));%
  % sprintf('mean');
  % disp(mean(ratio));%
  ratio_min_all(jfile) = min(ratio);
  ratio_max_all(jfile) = max(ratio);
  ratio_mean_all(jfile) = mean(ratio);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  roixy = roixy(masklocs,:);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  roixy = roixy(no_negs,:);
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  allxc = coords_file(:,1)-roixy(:,1)*pixelsize;
  allyc = coords_file(:,2)-roixy(:,2)*pixelsize;
  numlocsfile = size(coords_file,1);



  clear localizations mu

  % add data to PSF matrix
  runinserial = 1;
  if runinserial
    numworkersarg = Inf;
    p = gcp('nocreate');
    if isempty(p)
      myparpool = parpool;
    end
  else
    numworkersarg = 0;
  end

  parfor (ix = 1:Npatch,numworkersarg)    
    for iy = 1:Npatch
      fprintf('shift and sum ROI data patch (ix,iy)=(%i,%i).\n',ix,iy)
      for iz = 1:Nz
        % positionmask = (coords_file(:,1)>Xpatch(ix))&(coords_file(:,1)<Xpatch(ix+1))&...
        %                (coords_file(:,2)>Ypatch(iy))&(coords_file(:,2)<Ypatch(iy+1))&...
        %                (coords_file(:,3)>Zpatch(iz))&(coords_file(:,3)<Zpatch(iz+1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Include boundary points?
        positionmask = (coords_file(:,1)>=Xpatch(ix))&(coords_file(:,1)<=Xpatch(ix+1))&...
                       (coords_file(:,2)>=Ypatch(iy))&(coords_file(:,2)<=Ypatch(iy+1))&...
                       (coords_file(:,3)>=Zpatch(iz))&(coords_file(:,3)<=Zpatch(iz+1));
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        coords_select = [allxc,allyc];
        coords_select = coords_select(positionmask,:);
        coords_select = coords_select'/pixelsize/upsfac;
        allspots_select = allspots(:,:,positionmask);
        allmodelPSFs_select = allmodelPSFs(:,:,positionmask);
        allxc_select = allxc(positionmask);
        allyc_select = allyc(positionmask);

        allbg_select = bg(positionmask);
        allNph_select = Nph(positionmask);

        if isempty(coords_select)
          numlocsselect = 0;
        else
          numlocsselect = size(coords_select,2);
          bg_all(iz,ix,iy) = bg_all(iz,ix,iy)+sum(allbg_select);% store the total background count
          Nph_all(iz,ix,iy) = Nph_all(iz,ix,iy)+sum(allNph_select);% store the total photon count
          numlocmatrix(iz,ix,iy) = numlocmatrix(iz,ix,iy)+numlocsselect; % keep count of #ROIs per bin
          PSF_shift_ups = sum_shift_ups_PSFs(allspots_select,coords_select,upsfac); % shift ROI data to center and upsample
          PSFdatamatrix(:,:,iz,ix,iy) = PSFdatamatrix(:,:,iz,ix,iy)+PSF_shift_ups; % add ROI data
          PSF_shift_ups = sum_shift_ups_PSFs(allmodelPSFs_select,coords_select,upsfac); % shift model PSF data to center and upsample
          PSFmodelmatrix(:,:,iz,ix,iy) = PSFmodelmatrix(:,:,iz,ix,iy)+PSF_shift_ups; % add model PSF data
          allchisquares = (allspots_select-allmodelPSFs_select).^2./allmodelPSFs_select; % all chi-square values across the ROI
          chisquarematrix(:,:,iz,ix,iy) = chisquarematrix(:,:,iz,ix,iy)+sum(allchisquares,3);
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %% Here we calculate the scaled PSFs (with respect to intensity and background)
          % we need to apply 'PSF = (mu-bg)/Nph' over each individual
          % selection: so we reshape the bg and Nph 
          % this doesn't work properly for the datamatrix: maybe it doesn't
          % properly localize the points?
          % Norm_PSFs =(allspots_select - reshape(allbg_select, 1, 1, [])) ./ reshape(allNph_select, 1, 1, []);% from raw data to PSF (doesn't work ;_;)
          % % We cannot have negative values present in sum_shift_ups_PSFs:
          % % it messes with the centering. Either clip the negative values, or
          % % shift the minimum to 0 (or find another centering method).
          % min_PSFs = (min(min(Norm_PSFs, [], 1), [], 2))
          % Norm_PSFs = Norm_PSFs - min_PSFs;
          PSF_shift_ups = sum_shift_ups_PSFs_mod(allspots_select,coords_select,upsfac,allbg_select,allNph_select);
          % mu = PSF*Nph+bg, PSF = (mu-bg)/Nph (should also work for
          % allspots, mu and allspots semm to have the same size...
          PSF_Norm(:,:,iz,ix,iy) = PSF_Norm(:,:,iz,ix,iy)+PSF_shift_ups;
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
      end
    end
  end
end
sprintf('min: %i',min(ratio_min_all))
sprintf('max: %i',max(ratio_max_all))
sprintf('mean: %i',mean(ratio_mean_all))
%Normalize the PSF!
for ix = 1:Npatch
    for iy = 1:Npatch
        for iz = 1:Nz    
            PSF_Norm(:,:,iz,ix,iy) = PSF_Norm(:,:,iz,ix,iy)/numlocmatrix(iz,ix,iy);
        end
    end
end
% normalize chi-square values
for ix = 1:Npatch
  for iy = 1:Npatch
    for iz = 1:Nz
      chisquarematrix(:,:,iz,ix,iy) = chisquarematrix(:,:,iz,ix,iy)/numlocmatrix(iz,ix,iy);
    end
  end
end

%%
% save data

switch datacase
  case 'Lidke3D'
    savedatadir = 'C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\';
  case 'Fu3D'
    savedatadir = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu\data_analysis\analysis_vectorfit\PSF from data\';
end
%% I modified the name from 'PSFdata_'
filename = strcat(savedatadir,'PSFdatamatrix_no_negative_Nph_x7_',datacase);
save(filename,'PSFdatamatrix','PSFmodelmatrix','numlocmatrix','chisquarematrix','Xpatch','Ypatch','Zpatch',"PSF_Norm",'bg_all','Nph_all')

%%
% plot through focus PSF, data and model side by side

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

% create movie object
moviedir = 'C:\Users\Naam\Documents\BEP\BEP Code Roman\Data_BEP_Roman\All_files_data\';
writerObjPSF = VideoWriter(strcat(moviedir,'throughfocusPSF_matrix_clipped',datacase,'.avi'));
writerObjPSF.FrameRate = 1;
open(writerObjPSF);

for iz = 1:Nz
  fignum = iz;
  figure(fignum)
  figpos = [80 40 1440 720];
  set(gcf,'Position',figpos);
  tiledlayout(Npatch,2*Npatch,'TileSpacing','compact','Padding','compact')
  for ix = 1:Npatch
    for iy = 1:Npatch
      nexttile
      imagesc(PSFdatamatrix(:,:,iz,ix,iy))
      axis square
      axis off
      axis tight
      colormap parula
      nexttile
      imagesc(PSFmodelmatrix(:,:,iz,ix,iy))
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
  for ix = 1:Npatch
    for iy = 1:Npatch
      numlocstring = strcat('K=',num2str(numlocmatrix(iz,ix,iy)));
      locnumann{ix,iy} = annotation('textbox',[0.055+(iy-1)*0.229 0.91-0.23*(ix-1) 0.07 0.03],'String',numlocstring,'FontSize',12,'Edgecolor','none','Color','white');
    end
  end
  
  frame = getframe(gcf);
  writeVideo(writerObjPSF,frame);
end

close(writerObjPSF);

%%
% plot through focus PSF, bias between data and model

close all

% scale bar parameters
scalebarlength = 1;
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

locnumann = cell(Npatch,Npatch); %initialize cell array for annotation

% % position panels (x=rows=top-down, y=columns=left-right)
% panelposx = [0.71 0.49 0.27 0.05];
% panelposy = [0.03 0.26 0.49 0.73];
% panelsize = 130;

% create movie object
writerObjPSFbias = VideoWriter(strcat(moviedir,'throughfocusPSFbias_matrix_clipped',datacase,'.avi'));
writerObjPSFbias.FrameRate = 1;
open(writerObjPSFbias);

for iz = 1:Nz
  fignum = Nz+iz;
  figure(fignum)
  figpos = [80 40 1080 720];
  set(gcf,'Position',figpos);
  tiledlayout(Npatch,Npatch,'TileSpacing','compact','Padding','compact')
  for ix = 1:Npatch
    for iy = 1:Npatch
      nexttile
      imagesc(upsfac^2*(PSFdatamatrix(:,:,iz,ix,iy)-PSFmodelmatrix(:,:,iz,ix,iy))/numlocmatrix(iz,ix,iy))
      axis square
      axis off
      axis tight
      colormap parula
      colorbar
    end
  end
  sgtitle(strcat('z = ',num2str(zz(iz)),' nm'),'FontSize',20)
%   if iz==1
    annotation('rectangle',[0.11 0.75 (0.13/2.13)*scalebarlength 0.01],'FaceColor','white','Color','white');
%     annotation('rectangle',[0.10 0.75 0.13 0.01],'FaceColor','white','Color','white');
    annotation('textbox',[0.11 0.76 0.12 0.04],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
%   end
  for ix = 1:Npatch
    for iy = 1:Npatch
      annotation('textbox',[0.135+(iy-1)*0.234 0.948-0.23*(ix-1) 0.07 0.03],'String',strcat('(',num2str(ix),',',num2str(iy),')'),'FontSize',12,'Edgecolor','none','Color','black');
    end
  end
  
  frame = getframe(gcf);
  writeVideo(writerObjPSFbias,frame);
end

close(writerObjPSFbias);

%%

close all

switch datacase
  case 'Lidke3D'
    minlev = 0.8;
    maxlev = 2.0;
  case 'Fu3D'
    minlev = 0.8;
    maxlev = 2.8;
end

% create movie object for chi-square images
writerObjchisq = VideoWriter(strcat(moviedir,'chisquares_matrix_clipped',datacase,'.avi'));
writerObjchisq.FrameRate = 1;
open(writerObjchisq);

for iz = 1:Nz
  fignum = 2*Nz+iz;
  figure(fignum)
  figpos = [80 40 1080 720];
  set(gcf,'Position',figpos);
  tiledlayout(Npatch,Npatch,'TileSpacing','compact','Padding','compact')
  for ix = 1:Npatch
    for iy = 1:Npatch
      nexttile
      imagesc(chisquarematrix(:,:,iz,ix,iy))
      set(gca,'CLim',[minlev maxlev])
      axis square
      axis off
      axis tight
      colormap parula
      colorbar
    end
  end
  sgtitle(strcat('z = ',num2str(zz(iz)),' nm'),'FontSize',20)
%   if iz==1
    annotation('rectangle',[0.11 0.75 (0.13/2.13)*scalebarlength 0.01],'FaceColor','white','Color','white');
%     annotation('rectangle',[0.10 0.75 0.13 0.01],'FaceColor','white','Color','white');
    annotation('textbox',[0.11 0.76 0.12 0.04],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
%   end
  for ix = 1:Npatch
    for iy = 1:Npatch
      annotation('textbox',[0.145+(iy-1)*0.229 0.945-0.23*(ix-1) 0.07 0.03],'String',strcat('(',num2str(ix),',',num2str(iy),')'),'FontSize',12,'Edgecolor','none','Color','black');
    end
  end
    
  frame = getframe(gcf);
  writeVideo(writerObjchisq,frame);
end

close(writerObjchisq);

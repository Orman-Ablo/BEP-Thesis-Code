% This script is for linking localization data, in order to do a run length
% analysis and to make an assessment of localization precision.
%
% Sjoerd Stallinga, TUD, 2019-2023

close all
clear all

%%
% parameter settings, read-in data

fprintf('load unlinked localization data ...\n')

% read-in data
%datadir = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\localizations_after_drift_correction\';

% locdatafilename = 'localizations_27sept.mat'; 
% locdatafilename = 'localizations_const_astigmatism_27sept';
%locdatafilename = 'localizations_gain12_3okt';
%filename = strcat(datadir,locdatafilename);

%load(filename,'localizations','localizations_headers','params')

% parameters
pixelsize = params.pixelsize; % pixelsize in nm
minrunlength = 2; % minimum run length in the statistical evaluation
maxrunlength = 9; % maximum run length in the statistical evaluation
maxlocxy = 60; % max. for histogram
maxlocz = 100; % max. for histogram

% setup array for multiframeconnect
numlocs = size(localizations,1);
coords = localizations(:,2:4);
framenumber = localizations(:,5);
CRLBs = localizations(:,6:8);
Nphs = localizations(:,9);
bgs = localizations(:,10);
runlengths_in = 1*ones(numlocs,1);
locindex = (1:numlocs)';
% clear localizations
% loc = saveloc.loc;
% coords = cat(2,loc.xnm,loc.ynm,loc.znm);
% numlocs = size(coords,1);
% framenumber = loc.frame;
% CRLBs = cat(2,loc.locprecnm,loc.locprecnm,2*loc.locprecnm);
% Nphs = loc.phot;
% bgs = loc.bg;
% runlengths_in = 1*ones(numlocs,1);
% locindex = (1:numlocs)';

% filter for outlier fits
Nphmin = 100;
Nphmax = 1e4;
bgmin = -1;
bgmax = 150;
locfilter = (Nphs>Nphmin)&(Nphs<Nphmax)&(bgs>bgmin)&(bgs<bgmax)&...
            (CRLBs(:,1)<maxlocxy)&(CRLBs(:,2)<maxlocxy)&(CRLBs(:,3)<maxlocz);


%%
% analyze statistics of unlinked localizations
% use this to assess reasonable limits for CRLB, bg and Nph

% check statistics of initial localization events
scrsz = get(0,'ScreenSize');
figure('Position',[1*scrsz(4)/8 1*scrsz(4)/8 3*scrsz(3)/4 3*scrsz(4)/4]);
nbins = 130;
subplot(2,3,1);
histogram(CRLBs(locfilter,1),nbins);
title('CRLB x (nm)');
xlim([0 maxlocxy])
subplot(2,3,2);
histogram(CRLBs(locfilter,2),nbins);
title('CRLB y (nm)');
xlim([0 maxlocxy])
subplot(2,3,3);
histogram(CRLBs(locfilter,3),nbins);
title('CRLB z (nm)');
xlim([0 maxlocz])
subplot(2,3,4);
histogram(Nphs(locfilter),nbins);
title('signal photon count');
xlim([0 Nphmax])
subplot(2,3,5);
hist(bgs(locfilter),nbins);
title('background photons/pixel');
xlim([0 bgmax])

%%

fprintf('linking localization data ...\n')

% apply multiframe connect, localizations (xi,yi,ti) and (xj,yj,tj)
% (xy-coordinates and frame index) are merged if below conditions are all met:
% |xi-xj|<maxdist
% |yi-yj|<maxdist
% sqrt(xi-xj)^2+yi-yj)^2)<maxdist
% |xi-xj|<sigmafac*max(Delta_xi,Delta_xj)
% |yi-yj|<sigmafac*max(Delta_yi,Delta_yj)
% |ti-tj|<=nframes
%
sigmafac = 3; %std factor in determining the merger, sigmafac=3 is the recommended value, if it is too low only locs that are already very close are merged leading to an underestimation of the loc error, if it is too large locs are merged that are not from the same emitter
maxdistxy = 50; % max. distance factor (in same units as xy-coordinates, so nm at this point, set this to estimated half distance to nearest neighbor binding site)
maxdistz = 100; % max. distance factor (in same units as xy-coordinates, so nm at this point, set this to estimated half distance to nearest neighbor binding site)
nframes = 1; % max. number of consecutive frames that is checked, i.e. nframes-1 "dark" frames allowed
alldata = double([coords framenumber CRLBs.^2 CRLBs Nphs bgs runlengths_in coords coords.^2 locindex]);
%alldata = single([coords framenumber CRLBs.^2 CRLBs Nphs bgs runlengths_in coords coords.^2 locindex]);
alldata = alldata(locfilter,:);
alldata_connect = cMultiFrameConnect_allcalcs3D(alldata,sigmafac,maxdistxy,maxdistz,nframes);

% new variable names
numlocs_connect = size(alldata_connect,1);
coords_connect = alldata_connect(:,1:3);
index_connect = alldata_connect(:,4);
rstd_connect = sqrt(alldata_connect(:,5:7));
CRLBs_connect = alldata_connect(:,8:10);
Nphs_connect = alldata_connect(:,11);
bg_connect = alldata_connect(:,12);
runlengths = alldata_connect(:,13);
rsum = alldata_connect(:,14:16);
rsqsum = alldata_connect(:,17:19);

%clear alldata_connect

runlengthfilter = (runlengths>=minrunlength)&(runlengths<=maxrunlength);

% compute bias and localization error
rbias = zeros(numlocs_connect,3); % bias (= arithmetic mean-weighted average) and unbiased sample variance of x and y
rmse = zeros(numlocs_connect,3);
for kk = 1:3
  rbias(:,kk) = rsum(:,kk)./runlengths-coords_connect(:,kk);
  rmse(:,kk) = (rsqsum(:,kk)-(rsum(:,kk).^2)./runlengths)./(runlengths-1);
  rmse(runlengths==1,kk) = CRLBs_connect(runlengths==1,kk).^2;  % fix for unlinked runlength=1 events
end
locunc_connect = sqrt(abs(rmse));
for kk = 1:3
  locunc_connect(locunc_connect(:,kk)==0,kk) = mean(locunc_connect(runlengthfilter,kk)); % set loc.unc. for unmerged events equal to mean of loc. unc.
end

mean_locunc = mean(locunc_connect(runlengthfilter,:));
std_locunc = std(locunc_connect(runlengthfilter,:));
mean_CRLB = mean(CRLBs_connect(runlengthfilter,:));
std_CRLB = std(CRLBs_connect(runlengthfilter,:));

alllabelsxyz = {'x','y','z'};
for kk = 1:3
  fprintf(strcat('mean localization error',32,alllabelsxyz{kk},' = %5.2f nm\n'),mean_locunc(kk));
  fprintf(strcat('std localization error',32,alllabelsxyz{kk},' = %5.2f nm\n'),std_locunc(kk));
  fprintf(strcat('mean CRLB',32,alllabelsxyz{kk},' = %5.2f nm\n'),mean_CRLB(kk));
  fprintf(strcat('std CRLB',32,alllabelsxyz{kk},' = %5.2f nm\n'),std_CRLB(kk));
end

%%
figure
histogram(locunc_connect(runlengthfilter,1:2))
xlabel('localization precision x,y (nm)')
ax=gca;
ax.FontSize = 20;
%%
% check statistics of merged localization events

% show histograms
scrsz = get(0,'ScreenSize');
figure('Position',[1*scrsz(4)/8 1*scrsz(4)/8 3*scrsz(3)/4 3*scrsz(4)/4]);
nbins = 130;
allCRLBtitles = {'CRLB x (nm)','CRLB y (nm)','CRLB z (nm)'};
allxlims = [maxlocxy maxlocxy maxlocz];

for kk = 1:3
subplot(2,3,kk);
histogram(CRLBs_connect(runlengthfilter,kk),nbins);
title(allCRLBtitles{kk});
xlim([0 allxlims(kk)])
end
subplot(2,3,4);
histogram(Nphs_connect(runlengthfilter),nbins);
title('signal photon count');
xlim([0 Nphmax])
subplot(2,3,5);
histogram(bg_connect(runlengthfilter),nbins);
title('background photons/pixel');
xlim([0 4*bgmax])

%%
% compute on-time/average run length

% check run length distribution
alllengths = (1:maxrunlength)';
allnumlocs = zeros(maxrunlength,1);
for jj = 1:maxrunlength
   cinds = runlengths==jj;
   allnumlocs(jj) = sum(cinds);
end

% fit on-time
yy = allnumlocs(minrunlength:maxrunlength);
yy = yy/sum(yy);
yy = log(yy);
xx = (minrunlength:maxrunlength)';
fitorder = 1;
ccs = polyfit(xx,yy,fitorder);
ontime = -1/ccs(1);
fprintf('on-time fitted from exponential distribution = %3.2f\n',ontime);

% average run length
allnumlocs_norm = allnumlocs/sum(allnumlocs);
avrunlength = sum(alllengths.*allnumlocs_norm);
offsval = mean(allnumlocs_norm.*exp(alllengths/avrunlength));
fprintf('on-time averaged over all events = %3.2f\n',avrunlength);

% plot
figure
probrunlength = allnumlocs_norm(minrunlength:maxrunlength)/sum(allnumlocs_norm(minrunlength:maxrunlength));
basevaluepar = floor(2*min(log(probrunlength)/log(10)))/2;
barwd = 0.4;
bar(alllengths(minrunlength:maxrunlength),log(probrunlength)/log(10),barwd,'r','BaseValue',basevaluepar)
box on
hold on
plot(alllengths,(ccs(1)*alllengths+ccs(2))/log(10),'--k','LineWidth',2)
xlabel('on time [frames]','FontSize',12)
ylabel('log(probability)/log(10)','FontSize',12)
legend('data','exponential fit')
ylim([basevaluepar 0])

%%
% loop over regions in FOV to assess variations in performance

% define FOV patches
Npatch = 1;

minxyz = min(coords);
maxxyz = max(coords);
Xpatch = minxyz(1)+(0:Npatch)*(maxxyz(1)-minxyz(1))/Npatch;
Ypatch = minxyz(2)+(0:Npatch)*(maxxyz(2)-minxyz(2))/Npatch;

allmean_bg = zeros(Npatch,Npatch);
allstd_bg = zeros(Npatch,Npatch);
allmean_Nph = zeros(Npatch,Npatch);
allstd_Nph = zeros(Npatch,Npatch);
allmean_locunc = zeros(3,Npatch,Npatch);
allstd_locunc = zeros(3,Npatch,Npatch);
allmean_CRLB = zeros(3,Npatch,Npatch);
allstd_CRLB = zeros(3,Npatch,Npatch);
allmean_locunc_runl = zeros(3,maxrunlength,Npatch,Npatch);
allstd_locunc_runl = zeros(3,maxrunlength,Npatch,Npatch);
allmean_CRLB_runl = zeros(3,maxrunlength,Npatch,Npatch);
allstd_CRLB_runl = zeros(3,maxrunlength,Npatch,Npatch);

for ii = 1:Npatch
  for jj = 1:Npatch
    FOVfilter = (coords(:,1)>Xpatch(ii))&(coords(:,1)<Xpatch(ii+1))&...
                (coords(:,2)>Ypatch(jj))&(coords(:,2)<Ypatch(jj+1));
    FOVfilter_connect = (coords_connect(:,1)>Xpatch(ii))&(coords_connect(:,1)<Xpatch(ii+1))&...
                        (coords_connect(:,2)>Ypatch(jj))&(coords_connect(:,2)<Ypatch(jj+1));

    % bg and Nph across patches
    allmean_bg(ii,jj) = mean(bgs((bgs<bgmax)&FOVfilter));
    allstd_bg(ii,jj) = std(bgs((bgs<bgmax)&FOVfilter));
    allmean_Nph(ii,jj) = mean(Nphs((Nphs<Nphmax)&FOVfilter));
    allstd_Nph(ii,jj) = std(Nphs((Nphs<Nphmax)&FOVfilter));

    % make one-to-one correspondence of CRLBs of initial fits and measured
    % localization precision from the multiframeconnect, these histograms
    % should correspond with similar mean and median

    scrsz = get(0,'ScreenSize');
    figure('Position',[1*scrsz(4)/8 1*scrsz(4)/8 3*scrsz(3)/4 3*scrsz(4)/4]);
    nbins = 40;
    alllocunctitles = {'localization error in x (nm)','localization error in y (nm)','localization error in z (nm)'};

    for kk = 1:3
      subplot(2,3,kk);
      histogram(CRLBs((CRLBs(:,kk)<maxlocxy)&FOVfilter,kk),nbins);
      title(allCRLBtitles{kk});
      xlim([0 allxlims(kk)])
      subplot(2,3,3+kk);
      histogram(locunc_connect((locunc_connect(:,kk)<allxlims(kk))&FOVfilter_connect&runlengthfilter,kk),nbins);
      title(alllocunctitles{kk});
      xlim([0 allxlims(kk)])
    end
    
    % compute median, mean and std of CRLB and locunc
    for kk = 1:3
      locunc_median = median(locunc_connect(FOVfilter_connect&runlengthfilter&(locunc_connect(:,kk)<allxlims(kk)),kk));
      locunc_mean = mean(locunc_connect(FOVfilter_connect&runlengthfilter&(locunc_connect(:,kk)<allxlims(kk)),kk));
      locunc_std = std(locunc_connect(FOVfilter_connect&runlengthfilter&(locunc_connect(:,kk)<allxlims(kk)),kk));
      allmean_locunc(kk,ii,jj) = locunc_mean;
      allstd_locunc(kk,ii,jj) = locunc_std;
      
      fprintf(strcat('median rms localization error',32,alllabelsxyz{kk},' = %3.2f nm\n'),locunc_median);
      fprintf(strcat('mean rms localization error',32,alllabelsxyz{kk},' = %3.2f nm\n'),locunc_mean);
      fprintf(strcat('std rms localization error',32,alllabelsxyz{kk},' = %3.2f nm\n'),locunc_std);
    
      CRLB_median = median(CRLBs((CRLBs(:,kk)<maxlocxy)&FOVfilter,kk));
      CRLB_mean = mean(CRLBs((CRLBs(:,kk)<maxlocxy)&FOVfilter,kk));
      CRLB_std = std(CRLBs((CRLBs(:,kk)<maxlocxy)&FOVfilter,kk));
      allmean_CRLB(kk,ii,jj) = CRLB_mean;
      allstd_CRLB(kk,ii,jj) = CRLB_std;
      
      fprintf(strcat('median CRLB',32,alllabelsxyz{kk},' = %3.2f nm\n'),CRLB_median);
      fprintf(strcat('mean CRLB',32,alllabelsxyz{kk},' = %3.2f nm\n'),CRLB_mean);
      fprintf(strcat('std CRLB',32,alllabelsxyz{kk},' = %3.2f nm\n'),CRLB_std);
    end  
        
    % compute and plot ratio rms localization error to CRLB as a function of
    % run length, same plot as shown in Rieger/Stallinga, ChemPhysChem20014

    mean_locunc_runl = zeros(3,maxrunlength);
    mean_CRLB_runl = zeros(3,maxrunlength);
    mean_precisionratio = zeros(3,maxrunlength);
    std_locunc_runl = zeros(3,maxrunlength);
    std_CRLB_runl = zeros(3,maxrunlength);
    std_precisionratio = zeros(3,maxrunlength);
    for jr = minrunlength:maxrunlength
      for kk = 1:3
        CRLBxyz_tmp = CRLBs_connect(FOVfilter_connect&(runlengths==jr),kk);
        locuncxyz_tmp = locunc_connect(FOVfilter_connect&(runlengths==jr),kk);
        precisionratio_xyz = locuncxyz_tmp./CRLBxyz_tmp;
        mean_CRLB_runl(kk,jr) = mean(CRLBxyz_tmp);
        mean_locunc_runl(kk,jr) = mean(locuncxyz_tmp);
        mean_precisionratio(kk,jr) = mean(precisionratio_xyz);
        std_CRLB_runl(kk,jr) = std(CRLBxyz_tmp);
        std_locunc_runl(kk,jr) = std(locuncxyz_tmp);
        std_precisionratio(kk,jr) = std(precisionratio_xyz);
      end
    end

    % store data in array
    allmean_locunc_runl(:,:,ii,jj) = mean_locunc_runl;
    allstd_locunc_runl(:,:,ii,jj) = std_locunc_runl;
    allmean_CRLB_runl(:,:,ii,jj) = mean_CRLB_runl;
    allstd_CRLB_runl(:,:,ii,jj) = std_CRLB_runl;
    
    % plot overall xyz localization error as a function of run length
    scrsz = get(0,'ScreenSize');
    figure('Position',[1*scrsz(4)/8 1*scrsz(4)/8 3*scrsz(3)/4 3*scrsz(4)/4]);

    allylims = [30 30 50];
    allloclegs = {'loc. error x','loc. error y','loc. error z'};
    allCRLBlegs = {'CRLB x','CRLB y','CRLB z'};
    alllabelratios = {'rms loc. error / CRLB in x','rms loc. error / CRLB in y','rms loc. error / CRLB in z'};
    
    for kk = 1:3
      subplot(2,3,kk)
      % set(gcf,'units','pixels');
      % set(gcf,'Position',[100 100 450 300]);
      hold on
      box on
      errorbar((minrunlength:maxrunlength)-0.15,mean_locunc_runl(kk,minrunlength:maxrunlength),std_locunc_runl(kk,minrunlength:maxrunlength),'sb','MarkerSize',8,'LineWidth',2)
      errorbar((minrunlength:maxrunlength)+0.15,mean_CRLB_runl(kk,minrunlength:maxrunlength),std_CRLB_runl(kk,minrunlength:maxrunlength),'or','MarkerSize',8,'LineWidth',2)
      xlabel('on time [frames]')
      ylabel('rms error [nm]')
      legend(allloclegs{kk},allCRLBlegs{kk},'FontSize',14,'Box','off','NumColumns',2)
      xlim([minrunlength-1 maxrunlength+1])
      ylim([0 allylims(kk)])
      % xticks([(minrunlength-1):(maxrunlength+1)])
      set(gca,'FontSize',14)
      
      subplot(2,3,3+kk)
      hold on
      box on
      errorbar((minrunlength:maxrunlength),mean_precisionratio(kk,minrunlength:maxrunlength),std_precisionratio(kk,minrunlength:maxrunlength),'or','MarkerSize',8,'LineWidth',2)
      plot([minrunlength-1 maxrunlength+1],[1 1],'k--')
      xlabel('on time [frames]')
      ylabel(alllabelratios{kk})
      xlim([minrunlength-1 maxrunlength+1])
      ylim([0 2])
      set(gca,'FontSize',14)
    end
    
  end
end

%%

% plot found localization errors across FOV
scrsz = get(0,'ScreenSize');
figure('Position',[1*scrsz(4)/8 1*scrsz(4)/8 3*scrsz(3)/4 3*scrsz(4)/4]);

climsxyz = [8 16; 8 16;16 24];

for kk = 1:3
  subplot(2,3,kk)
  imagesc(squeeze(allmean_CRLB(kk,:,:)),climsxyz(kk,:))
  axis square
  title(allCRLBtitles{kk}) 
  colorbar
  
  subplot(2,3,3+kk)
  imagesc(squeeze(allmean_locunc(kk,:,:)),climsxyz(kk,:))
  axis square
  title(alllocunctitles{kk}) 
  colorbar
end

% plot for bg and Nph across FOV
figure

subplot(1,2,1)
imagesc(allmean_bg)
axis square
title('bg (photons/pix)') 
colorbar

subplot(1,2,2)
imagesc(allmean_Nph)
axis square
title('Nph (photons)') 
colorbar

%%
% save data linked localization

fprintf('save linked localization data ...\n')

datadir = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Fu_data_analysis\localizations_after_linking\';
filename = strcat(datadir,locdatafilename);

save(filename,'coords_connect','index_connect',...
              'locunc_connect','CRLBs_connect',...
              'Nphs_connect','bg_connect','runlengths',...
              'allmean_locunc','allstd_locunc','allmean_CRLB','allstd_CRLB',...
              'allmean_Nph','allmean_bg');
            
%%
% make SR-images of original unlinked and with linked localizations
% using histogram visualization

% make overview image, widefield equivalent
SRpixelsize = 100.0;
SRpixelsizez = 1000.0;

% unlinked localizations
xall_orig = coords(locfilter,1);
yall_orig = coords(locfilter,2);
zall_orig = coords(locfilter,3);
minxyz = min(coords);
maxxyz = max(coords);
minxyz(3) = -500;
maxxyz(3) = 500;
rangexyz = maxxyz-minxyz;
Nx = ceil(rangexyz(1)/SRpixelsize);
Ny = ceil(rangexyz(2)/SRpixelsize);
Nz = ceil(rangexyz(3)/SRpixelsizez);
Xedges = minxyz(1)+(0:Nx)*SRpixelsize;
Yedges = minxyz(2)+(0:Ny)*SRpixelsize;
Zedges = minxyz(3)+(0:Nz)*SRpixelsizez;

% get histogam
[im_hist,~,~,~,~] = histcounts2(xall_orig,yall_orig,Xedges,Yedges);

% 1 SR pixel Gaussian blurring
sigfac = 1.0;
im_filt = imgaussfilt(im_hist,sigfac);

% create figure
scrsz = get(0,'ScreenSize');
figure('Position',[1*scrsz(4)/8 1*scrsz(4)/8 3*scrsz(4)/4 3*scrsz(4)/4]);
xrange = [minxyz(1) maxxyz(1)];
yrange = [minxyz(2) maxxyz(2)];
CLim = [0 100];
imagesc(xrange,yrange,im_filt,CLim)
axis square
colorbar
colormap hot
xlabel('y [nm]')
ylabel('x [nm]')
set(gca,'XAxisLocation','top')

% create crop for swift inspection
xcrop = [1e4 5e4];
ycrop = [7e4 11e4];
% xcrop = [13.5e4 16.5e4]; % cell bottom right
% ycrop = [13.5e4 16.5e4];
SRpixelsize = 5.0;
SRpixelsizez = 10.0;
Nx = ceil((xcrop(2)-xcrop(1))/SRpixelsize);
Ny = ceil((ycrop(2)-ycrop(1))/SRpixelsize);
Xedges = xcrop(1)+(0:Nx)*SRpixelsize;
Yedges = ycrop(1)+(0:Ny)*SRpixelsize;

% get histogam
cropmask = (xall_orig>xcrop(1))&(xall_orig<xcrop(2))&(yall_orig>ycrop(1))&(yall_orig<ycrop(2));
[im_hist,~,~,~,~] = histcounts2(xall_orig(cropmask),yall_orig(cropmask),Xedges,Yedges);

% 1 SR pixel Gaussian blurring
sigfac = 1.0;
im_filt = imgaussfilt(im_hist,sigfac);

% create figure
scrsz = get(0,'ScreenSize');
figure('Position',[1*scrsz(4)/8 1*scrsz(4)/8 3*scrsz(4)/4 3*scrsz(4)/4]);
xrange = [xcrop(1) xcrop(2)];
yrange = [ycrop(1) ycrop(2)];
CLim = [0 1];
imagesc(xrange,yrange,im_filt,CLim)
set(gca,'XAxisLocation','top')
axis square
colorbar
colormap hot
xlabel('y [nm]')
ylabel('x [nm]')

% get histogam linked localizations
xall_connect = coords_connect(:,1);
yall_connect = coords_connect(:,2);
cropmask = (xall_connect>xcrop(1))&(xall_connect<xcrop(2))&(yall_connect>ycrop(1))&(yall_connect<ycrop(2));
[im_hist,~,~,~,~] = histcounts2(xall_connect(cropmask),yall_connect(cropmask),Xedges,Yedges);

% 1 SR pixel Gaussian blurring
sigfac = 1.0;
im_filt = imgaussfilt(im_hist,sigfac);

% create figure
scrsz = get(0,'ScreenSize');
figure('Position',[1*scrsz(4)/8 1*scrsz(4)/8 3*scrsz(4)/4 3*scrsz(4)/4]);
xrange = [xcrop(1) xcrop(2)];
yrange = [ycrop(1) ycrop(2)];
CLim = [0 1];
imagesc(xrange,yrange,im_filt,CLim)
set(gca,'XAxisLocation','top')
axis square
colorbar
colormap hot
xlabel('y [nm]')
ylabel('x [nm]')

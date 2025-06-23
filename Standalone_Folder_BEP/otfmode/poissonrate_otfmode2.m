function [mu,dmudtheta] = poissonrate_otfmode2(params,theta,OTF,intensity)
% Returns the Poisson-rates for all pixels and all first order derivatives
% w.r.t. the parameters theta.


% Modified version of 'poissonrate', to use the OTF for fitting 

K = params.K;
Mx = params.Mx;
My = params.My;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add in the 'compders' parameter for the OTF fit!
compders = params.compders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fitmodel = params.fitmodel;
numparams = params.numparams; % nr of local params

switch fitmodel
    case {'xy','xy-gamma'}
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
        
    case {'xyz','xyz-gamma'}
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
   
    case 'xy-azim'
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
        azim = theta(5);
        params.azim = azim;
        
    case 'xy-azim-pola'
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
        azim = theta(5);
        pola = theta(6);
        params.azim = azim;
        params.pola = pola;
        
    case 'xyz-azim-pola'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        azim = theta(6);
        pola = theta(7);
        params.azim = azim;
        params.pola = pola;
        
    case 'xy-azim-diffusion'
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
        azim = theta(5);
        params.azim = azim;
        params.g2 = theta(6);
        
    case 'xy-azim-pola-diffusion'
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
        azim = theta(5);
        pola = theta(6);
        params.azim = azim;
        params.pola = pola;
        params.g2 = theta(7);
        
    case 'xyz-azim-pola-diffusion'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        azim = theta(6);
        pola = theta(7);
        params.azim = azim;
        params.pola = pola;
        params.g2 = theta(8);
        
    case 'xyz-aberrations'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        params.aberrations(:,3) = theta(6:end);
    
    case 'xyz-tilt-aberrations'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        params.tilt = theta(6:7);
        params.aberrations(:,3) = theta(8:end);     
        
    case 'xyz-azim-pola-aberrations'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        azim = theta(6);
        pola = theta(7);
        params.azim = azim;
        params.pola = pola;
        params.aberrations(:,3) = theta(8:end);
        
    case 'xyz-azim-pola-diffusion-aberrations'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        azim = theta(6);
        pola = theta(7);
        params.azim = azim;
        params.pola = pola;
        params.g2 = theta(8);
        params.aberrations(:,3) = theta(9:end);
    otherwise
        error('fitmodel is not correct or not implemented.')
end

switch params.excitation
    case 'constant'
        P = 1;
        dPdazim = 0;
        dPdpola = 0;
    case 'zstack'
        P = ones(1,K);
        dPdazim = zeros(1,K);
        dPdpola = zeros(1,K);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the PSFs using the OTF
[PSF,PSFder] = get_psfs_derivatives_otfmode2(OTF,params,compders);
% disp(sum(PSF,'all'))

%% this produces the post-computation PSF and derivative scaling
pos = params.pos; % pixel position
Zpatch =params.Zpatch;
dz = Zpatch(2)-Zpatch(1); % z- interval size
z_emit = params.zemit/dz; % scaled pixel position
val = interp1(pos, intensity, z_emit, 'linear', 'extrap'); % linear interpolation for less computation power
scale = val/sum(PSF,'all');
PSF = PSF*scale;
% disp(sum(PSF,'all'))
PSFder(:,:,1,1:3) = PSFder(:,:,1,1:3)*scale; % scale the position derivatives

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get Poisson rate and derivatives
mu = zeros(Mx,My,K);
dmudtheta = zeros(Mx,My,K,numparams);

if strcmp(params.excitation,'zstack')

    if strcmp(fitmodel,'xy')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1:2) = Nph*PSFder(:,:,:,1:2);
        dmudtheta(:,:,:,3) = PSF;
        dmudtheta(:,:,:,4) = 1;
    end

    if strcmp(fitmodel,'xyz')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1:3) = Nph*PSFder(:,:,:,1:3);
        dmudtheta(:,:,:,4) = PSF;
        dmudtheta(:,:,:,5) = 1;
    end

    if strcmp(fitmodel,'xy-gamma')
        error('fitmodel xy-gamma not implemented for params.excitation=z-stack')
    end
    
    if strcmp(fitmodel,'xyz-aberrations')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1:3) = Nph*PSFder(:,:,:,1:3);
        dmudtheta(:,:,:,4) = PSF;
        dmudtheta(:,:,:,5) = 1;
        dmudtheta(:,:,:,6:end) = Nph*PSFder(:,:,:,4:end);
    end

    if strcmp(fitmodel,'xyz-tilt-aberrations')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1:3) = Nph*PSFder(:,:,:,1:3);
        dmudtheta(:,:,:,4) = PSF;
        dmudtheta(:,:,:,5) = 1;
        dmudtheta(:,:,:,6) = Nph*PSFder(:,:,:,4);
        dmudtheta(:,:,:,7) = Nph*PSFder(:,:,:,5);
        dmudtheta(:,:,:,8:end) = Nph*PSFder(:,:,:,6:end);
    end
    
    if strcmp(fitmodel,'xyz-azim-pola-aberrations')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1) = Nph*PSFder(:,:,:,1);
        dmudtheta(:,:,:,2) = Nph*PSFder(:,:,:,2);
        dmudtheta(:,:,:,3) = Nph*PSFder(:,:,:,3);
        dmudtheta(:,:,:,4) = PSF;
        dmudtheta(:,:,:,5) = 1/K;
        dmudtheta(:,:,:,6) = Nph*PSFder(:,:,:,4);
        dmudtheta(:,:,:,7) = Nph*PSFder(:,:,:,5);
        dmudtheta(:,:,:,8:end) = Nph*PSFder(:,:,:,6:end);
    end
    
    if strcmp(fitmodel,'xyz-azim-pola-diffusion-aberrations')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1) = Nph*PSFder(:,:,:,1);
        dmudtheta(:,:,:,2) = Nph*PSFder(:,:,:,2);
        dmudtheta(:,:,:,3) = Nph*PSFder(:,:,:,3);
        dmudtheta(:,:,:,4) = PSF;
        dmudtheta(:,:,:,5) = 1/K;
        dmudtheta(:,:,:,6) = Nph*PSFder(:,:,:,4);
        dmudtheta(:,:,:,7) = Nph*PSFder(:,:,:,5);
        dmudtheta(:,:,:,8) = Nph*PSFder(:,:,:,6);
        dmudtheta(:,:,:,9:end) = Nph*PSFder(:,:,:,7:end);
    end
    
    if strcmp(fitmodel,'xyz-azim-pola-diffusion')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1) = Nph*PSFder(:,:,:,1);
        dmudtheta(:,:,:,2) = Nph*PSFder(:,:,:,2);
        dmudtheta(:,:,:,3) = Nph*PSFder(:,:,:,3);
        dmudtheta(:,:,:,4) = PSF;
        dmudtheta(:,:,:,5) = 1/K;
        dmudtheta(:,:,:,6) = Nph*PSFder(:,:,:,4);
        dmudtheta(:,:,:,7) = Nph*PSFder(:,:,:,5);
        dmudtheta(:,:,:,8) = Nph*PSFder(:,:,:,6);
    end
    
else
    for k = 1:K
        % PSF model
        mu(:,:,k) = Nph*P(k)*PSF+Nbg/K;
        
        % get derivatives of Poisson rate w.r.t. fit parameters
        if strcmp(fitmodel,'xy-azim-pola') && K>1
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = P(k)*PSF;
            dmudtheta(:,:,k,4) = 1/K;
            dmudtheta(:,:,k,5) = Nph*P(k)*PSFder(:,:,3)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4)+Nph*dPdpola(k)*PSF;
        end
        
        if strcmp(fitmodel,'xy') || strcmp(fitmodel,'xy-gamma')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = P(k)*PSF;
            dmudtheta(:,:,k,4) = 1/K;
        end
        
        if strcmp(fitmodel,'xyz') || strcmp(fitmodel,'xyz-gamma')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = Nph*P(k)*PSFder(:,:,3);
            dmudtheta(:,:,k,4) = P(k)*PSF;
            dmudtheta(:,:,k,5) = 1/K;
        end

        if strcmp(fitmodel,'xy-azim')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = P(k)*PSF;
            dmudtheta(:,:,k,4) = 1/K;
            dmudtheta(:,:,k,5) = Nph*dPdazim(k)*PSF+Nph*P(k)*PSFder(:,:,3);
        end
        
        if strcmp(fitmodel,'xyz-azim-pola')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = Nph*P(k)*PSFder(:,:,3);
            dmudtheta(:,:,k,4) = P(k)*PSF;
            dmudtheta(:,:,k,5) = 1/K;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,7) = Nph*P(k)*PSFder(:,:,5)+Nph*dPdpola(k)*PSF;
        end
        
        if strcmp(fitmodel,'xy-azim-diffusion')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = P(k)*PSF;
            dmudtheta(:,:,k,4) = 1/K;
            dmudtheta(:,:,k,5) = Nph*P(k)*PSFder(:,:,3)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4);
        end
        
        if strcmp(fitmodel,'xy-azim-pola-diffusion')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = P(k)*PSF;
            dmudtheta(:,:,k,4) = 1/K;
            dmudtheta(:,:,k,5) = Nph*P(k)*PSFder(:,:,3)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4)+Nph*dPdpola(k)*PSF;
            dmudtheta(:,:,k,7) = Nph*P(k)*PSFder(:,:,5);
        end
        
        if strcmp(fitmodel,'xyz-azim-pola-diffusion')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = Nph*P(k)*PSFder(:,:,3);
            dmudtheta(:,:,k,4) = P(k)*PSF;
            dmudtheta(:,:,k,5) = 1/K;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,7) = Nph*P(k)*PSFder(:,:,5)+Nph*dPdpola(k)*PSF;
            dmudtheta(:,:,k,8) = Nph*P(k)*PSFder(:,:,6);
        end
        
    end
end
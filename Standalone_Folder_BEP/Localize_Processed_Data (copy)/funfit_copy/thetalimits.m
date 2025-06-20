% calculate theta estimation limits (max/min values)
function [thetamin,thetamax] = thetalimits(params,theta)

fitmodel = params.fitmodel;
zernikecoefsmax = 0.5*params.lambda*ones(1,size(params.aberrations,1));

Mx = params.Mx;
My = params.My;
roisizex = params.Mx*params.pixelsize;
roisizey = params.My*params.pixelsize;
xmin = -roisizex/2;
xmax = roisizex/2;
ymin = -roisizey/2;
ymax = roisizey/2;
zmin = params.zspread(1);
zmax = params.zspread(2);
azimmax = 2*pi;
polamax = pi;

if strcmp(fitmodel,'xy')
    thetamin = [xmin,ymin,theta(3)/10,theta(4)/10];
    thetamax = [xmax,ymax,2*theta(3),max(theta(3)/Mx/My/2,2*theta(4))];
elseif strcmp(fitmodel,'xyz') || strcmp(fitmodel,'xyz-gamma')
    thetamin = [xmin,ymin,zmin,theta(4)/10,theta(5)/10];
    thetamax = [xmax,ymax,zmax,2*theta(4),max(theta(4)/Mx/My/2,2*theta(5))];
elseif strcmp(fitmodel,'xy-azim')
    thetamin = [xmin,ymin,theta(3)/10,0,0];
    thetamax = [xmax,ymax,roisizey/2,2*theta(3),theta(3),azimmax];
elseif strcmp(fitmodel,'xy-azim-pola')
    thetamin = [xmin,ymin,theta(3)/10,0,0,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3),azimmax,polamax];
elseif strcmp(fitmodel,'xy-azim-pola-diffusion')
    thetamin = [xmin,ymin,theta(3)/10,0,0,-0,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3),azimmax,polamax,1];
elseif strcmp(fitmodel,'xyz-azim-pola')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax];
elseif strcmp(fitmodel,'xy-azim-diffusion')
    thetamin = [xmin,ymin,theta(3)/10,0,0,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3),azimmax,1];
elseif strcmp(fitmodel,'xyz-azim-pola-diffusion')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0,0];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax,1];
elseif strcmp(fitmodel,'xyz-aberrations')
    thetamin = [xmin,ymin,zmin,theta(4)/10,theta(5)/10,-zernikecoefsmax];
    thetamax = [xmax,ymax,zmax,2*theta(3),max(theta(3)/Mx/My/2,2*theta(4)),zernikecoefsmax];
elseif strcmp(fitmodel,'xyz-tilt-aberrations')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,-Inf,-Inf,-zernikecoefsmax];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),Inf,Inf,zernikecoefsmax];
elseif strcmp(fitmodel,'xyz-azim-pola-aberrations')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0,-zernikecoefsmax];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax,zernikecoefsmax];
elseif strcmp(fitmodel,'xyz-azim-pola-diffusion-aberrations')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0,0,-zernikecoefsmax];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax,1,zernikecoefsmax];
else 
    error('Fitmodel not correct or not implemented')
end

end

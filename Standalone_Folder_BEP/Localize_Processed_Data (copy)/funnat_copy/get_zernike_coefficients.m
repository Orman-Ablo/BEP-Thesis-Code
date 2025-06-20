function [zernikeCoefficients] = get_zernike_coefficients(xn, yn, gammas, params)
% This function calculates the Zernike coefficients, given gamma
% coefficients and x,y locations.
%
% INPUT
% - xn,yn: coordinates or meshgrid of points for which we want to determine
% the Zernike coefficients. Units: physical coordinates of the total FOV
% (um) where the origin is at the outer edge of pixel 1,1.
% - gamma: the gamma coefficients that determine how the Zernike
% Coefficients depend on x,y. 
%
% OUTPUT
% - zernikeCoefficients: an array or surface containing the zernike
% coeffcients in nm for each point (xn,yn).
%
% Autor: Isabel Droste, TU Delft, 2022, 2023

if size(xn,1) == 1
    xn = xn';
    yn = yn';
end
    
NATgammas = params.NATgammas;
orders2D = params.orders2D;
nr_of_zernike_modes = size(params.aberrations,1);
numgammas = length(NATgammas);

if size(xn,1) == 1
    xn = xn';
    yn = yn';
end

size_grid = size(xn);
zernikeCoefficients = zeros(size_grid(1),size_grid(2),nr_of_zernike_modes);

alllegendres2D = get_legendrepolynomials2D(orders2D,xn,yn);
legendre_normfac = sqrt((1+2*orders2D(:,1)).*(1+2*orders2D(:,2))); 

% loop over NATgammas, add up per zernike
jv = 1;
for jgam = 1:numgammas
    RR = NATgammas{jgam};
    multiples = length(RR)-1;
    for jm = 1:multiples
        NATvec = RR{jm};
        nn = NATvec(3);
        mm = NATvec(4);
        weight = NATvec(5);
        jzer = intersect(find(params.aberrations(:,1)==nn),find(params.aberrations(:,2)==mm)); % find the right Zernike in the list
        zernikeCoefficients(:,:,jzer) = zernikeCoefficients(:,:,jzer) + gammas(jgam)*weight*legendre_normfac(jv)*alllegendres2D(:,:,jv);
        jv = jv+1;
    end
    
end

if size_grid(2) == 1
    zernikeCoefficients = permute(zernikeCoefficients,[1 3 2]);
end
     
end
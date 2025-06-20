% This script can fit the gamma coefficients from Zernike coefficients
% fitted from bead data. It works for defocus, astigmatism, coma and
% spherical aberration.
%
% Author: Isabel Droste, TU Delft, 2024

function gammas = invertNAT(xn,yn,allzernikes,params)
Ncfg = numel(xn);

% Calculate the 2D Legendre polynomials
legendre_orders =  [ ...
    0,  0;
    1,  0;
    0,  1;
    1,  1;
    2,  0;
    0,  2;];
alllegendres2D = squeeze(get_legendrepolynomials2D(legendre_orders,xn,yn));
legendre_normfactor = sqrt((1+2*legendre_orders(:,1)).*(1+2*legendre_orders(:,2)))';
alllegendres2D = repmat(legendre_normfactor,[Ncfg,1]).*alllegendres2D;

% Build the Legendre matrices (Anm = L*gammas)
L_defocus = [alllegendres2D(:,1) alllegendres2D(:,2) alllegendres2D(:,3) alllegendres2D(:,5)+alllegendres2D(:,6)];
L_astigmatism = [zeros(Ncfg,1) alllegendres2D(:,1) alllegendres2D(:,3) alllegendres2D(:,2) alllegendres2D(:,4);...
                alllegendres2D(:,1) zeros(Ncfg,1) alllegendres2D(:,2) -alllegendres2D(:,3) (3/sqrt(5))*(alllegendres2D(:,5)-alllegendres2D(:,6));];
L_coma = [zeros(Ncfg,1) alllegendres2D(:,1) alllegendres2D(:,3); ...
            alllegendres2D(:,1) zeros(Ncfg,1) alllegendres2D(:,2);];
L_spherical = [alllegendres2D(:,1)];

% Build the Zernike vectors Anm
zernikes_defocus = allzernikes(intersect(find(params.aberrations(:,1)==2),find(params.aberrations(:,2)==0)),:)';
zernikes_astigmatism = [allzernikes(intersect(find(params.aberrations(:,1)==2),find(params.aberrations(:,2)==-2)),:) ...
                        allzernikes(intersect(find(params.aberrations(:,1)==2),find(params.aberrations(:,2)==2)),:)]';
zernikes_coma = [allzernikes(intersect(find(params.aberrations(:,1)==3),find(params.aberrations(:,2)==-1)),:) ...
                        allzernikes(intersect(find(params.aberrations(:,1)==3),find(params.aberrations(:,2)==1)),:)]';
zernikes_spherical = allzernikes(intersect(find(params.aberrations(:,1)==4),find(params.aberrations(:,2)==0)),:)';


% Find the least squares solution for gammas using the pseudoinverse
gammas_defocus = pinv(L_defocus)*zernikes_defocus;
gammas_astigmatism = pinv(L_astigmatism)*zernikes_astigmatism;
gammas_coma = pinv(L_coma)*zernikes_coma;
gammas_spherical = pinv(L_spherical)*zernikes_spherical;

gammas = [gammas_defocus; gammas_astigmatism; gammas_coma; gammas_spherical;];

end
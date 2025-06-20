function [PupilFunctionDerivative] = get_pupil_function_derivatives_nat(jgamma,zernike_normfac,allzernikes,PupilFunction,FOV_x,FOV_y,params)
%get_pupil_function_derivatives Calculate the pupil function derivatives
%with respect to the gamma variables.

% INPUT 
% jgamma = the number of the gamma variable from 1 until 14. 
% zernike_normfac = the norm factors of the zernike coefficients

% OUTPUT
% pupilfunctionDerivative = the derivative of the pupil function with
% respect to gamma variable jgamma


switch jgamma

    % Defocus: A(2,0)
    case 1
        PupilFunctionDerivative = (2*pi*1i*zernike_normfac(1)*allzernikes(:,:,1)/params.lambda).*PupilFunction;
    case 2
        PupilFunctionDerivative = (2*pi*1i*FOV_x*zernike_normfac(1)*allzernikes(:,:,1)/params.lambda).*PupilFunction;
    case 3
        PupilFunctionDerivative = (2*pi*1i*FOV_y*zernike_normfac(1)*allzernikes(:,:,1)/params.lambda).*PupilFunction;
    case 4
        PupilFunctionDerivative = (2*pi*1i*(FOV_x^2 + FOV_y^2)*zernike_normfac(1)*allzernikes(:,:,1)/params.lambda).*PupilFunction;
    
    % Astigmatism: A(2,-2) and A(2,2)
    case 5
        PupilFunctionDerivative = (2*pi*1i*zernike_normfac(2)*allzernikes(:,:,2)/params.lambda).*PupilFunction;
    case 6
        PupilFunctionDerivative = (2*pi*1i*zernike_normfac(3)*allzernikes(:,:,3)/params.lambda).*PupilFunction;
    case 7
        PupilFunctionDerivative = (2*pi*1i*(FOV_x*zernike_normfac(2)*allzernikes(:,:,2) - FOV_y*zernike_normfac(3)*allzernikes(:,:,3))/params.lambda).*PupilFunction;
    case 8
        PupilFunctionDerivative = (2*pi*1i*(FOV_y*zernike_normfac(2)*allzernikes(:,:,2) + FOV_x*zernike_normfac(3)*allzernikes(:,:,3))/params.lambda).*PupilFunction;
    case 9
        PupilFunctionDerivative = (2*pi*1i*(2*FOV_x*FOV_y*zernike_normfac(2)*allzernikes(:,:,2) + (FOV_x^2 - FOV_y^2)*zernike_normfac(3)*allzernikes(:,:,3))/params.lambda).*PupilFunction;

    % Coma: A(3,-1) and A(3,1)
    case 10
        PupilFunctionDerivative = (2*pi*1i*zernike_normfac(4)*allzernikes(:,:,4)/params.lambda).*PupilFunction;
    case 11
        PupilFunctionDerivative = (2*pi*1i*zernike_normfac(5)*allzernikes(:,:,5)/params.lambda).*PupilFunction;
    case 12
        PupilFunctionDerivative = (2*pi*1i*(FOV_y*zernike_normfac(4)*allzernikes(:,:,4) + FOV_x*zernike_normfac(5)*allzernikes(:,:,5))/params.lambda).*PupilFunction;

    % Spherical aberration: A(4,0)
    case 13
        PupilFunctionDerivative = (2*pi*1i*zernike_normfac(6)*allzernikes(:,:,6)/params.lambda).*PupilFunction;
        
    otherwise
        error('jgamma is not correct or not implemented.')
end 



end
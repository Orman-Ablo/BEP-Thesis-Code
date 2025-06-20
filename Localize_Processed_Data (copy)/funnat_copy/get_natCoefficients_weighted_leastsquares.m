function [RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6] = get_natCoefficients_weighted_leastsquares(xn,yn,zers,params)
% This function fits a Zernike surface using weighted least squares fitting
% where each data point is weighted with the CRLB value of the specific
% aberration for this data point. In this way, aberration that are determined
% more precisely, have a larger weight in the fitting. When using an
% ordinary least squares fitting, each data point has the same weight.

% Input
% xn, yn: the fitted x,y coordinates. This can be found in theta(1:2,:) as
% returned by function localization

% zers: the fitted zernike coefficients. This can be found in theta(6:end,:) as
% returned by function localization

% params: the fitting parameters that have been used as input for the
% localization function.


% Output
% RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6
% Each of these contain the coefficients that determine the low order
% polynomial that describes how the aberration coefficient varies over the
% field of view.



end
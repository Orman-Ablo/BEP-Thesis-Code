function [system_aberration_value] = get_system_aberration_value(gammas,params)
%get_system_aberration_value Calculates the system aberration value in
%mlambda from gammas and the multiplex factors.

multiplex_factors = params.multiplex_factors;
lambda = params.lambda;
system_aberration_value = sqrt(sum(multiplex_factors.*(gammas.^2)))*1e3/lambda;

end
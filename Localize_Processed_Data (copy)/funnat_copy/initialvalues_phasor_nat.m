function [thetainit] = initialvalues_phasor_nat(params)
% The initial values for theta global.
% Random values with a specific system aberration level 
%
%
% INPUT
% -params
%
% OUTPUT
% thetainit = struct with fields:
% thetainit.global
%
% Autor: Isabel Droste, TU Delft, 2022, 2023

%gammas_init = params.gammas_init;
numgammas = params.numgammas;
gammas_init_fixed = params.gammas_init;

% Initialize gammas randomly with a given system aberration level
gammas_init = (rand(numgammas,1)-0.5).*params.gammas_fitted;
scale_factor = params.aberration_level_init/get_system_aberration_value(gammas_init,params);
gammas_init = scale_factor*gammas_init;

%%%% test
%gammas_init = zeros(numgammas,1);
%%%%%%

% Replace values by pre-chosen initial value (astigmatism)
for i=find(gammas_init_fixed)
    gammas_init(i) = gammas_init_fixed(i);
end
% Make sure that gamma(5) is positive because
% want this solution and not the mirrored solution
gammas_init(5,:) = abs(gammas_init(5,:));

thetainit = gammas_init;

end

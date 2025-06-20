% calculate theta estimation limits (max/min values)
function [thetamin,thetamax] = thetalimits_nat(params)

%thetamin = 100*params.gammas_range(:,1);
%thetamax = 100*params.gammas_range(:,2);
numgammas = params.numgammas;
thetamin = -Inf*ones(numgammas,1); % select a better range
thetamax = Inf*ones(numgammas,1);

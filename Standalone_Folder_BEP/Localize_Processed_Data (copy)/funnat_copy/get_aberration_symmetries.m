function theta_global_symmetries = get_aberration_symmetries(theta_global,params)
%get_aberration_symmetries Calculates the gammas for the aberation
%symmetries
%   Calculates the gammas that lead to the same image as the given
%   theta_global. 
%   The gammas that belong to defocus, astigmatism and spherical aberration
%   are flipped, leading to at most 8 different sets of gammas.
%   Only the gammas that are fitted (params.gammas_fitted == true), 
%   are flipped sign.
gammas_fitted = params.gammas_fitted;

% Check which gammas are fitted. If at least one gamma that belongs to an
% aberration that is fitted, this aberration should be flipped.
flip_aberrations = [0 0 0];
fitted_gammas = find(gammas_fitted);

for k=1:size(fitted_gammas)
    i = fitted_gammas(k);
    switch num2str(i)
        case {'1', '2', '3', '4'}
            flip_aberrations(1) = 1;
        case {'5' '6' '7' '8' '9'}
            flip_aberrations(2) = 1;
        case{'14'}
            flip_aberrations(3) = 1;
        otherwise
            
    end          
end

% Create array of 0 (=no flip) and 1(=flip) for all possibilities (max 8).
flip_aberration_indices = find(flip_aberrations);
all_combs = dec2bin(0:2^(numel(flip_aberration_indices))-1)-'0';

% insert columns of zeros for the aberration that are not flipped.
for i=find(flip_aberrations==0)
    all_combs = [all_combs(:,1:i-1) zeros(size(all_combs,1),1) all_combs(:,i:end)];
end

% Create flipped solutions.
theta_global_symmetries = zeros([size(theta_global),size(all_combs,1)-1]);

for i=2:size(all_combs,1)
    gammas = theta_global;
    comb = all_combs(i,:);
    if comb(1)
        gammas(1:4) = -1*theta_global(1:4);
    end
    if comb(2)
        gammas(5:9) = -1*theta_global(5:9);
    end
    if comb(3)
        gammas(14) = -1*theta_global(14);
    end

    theta_global_symmetries(:,:,i-1) = gammas;

end

% todo :uitzondering. als beide astigmatismen alleen met een constante worden
% gefit, kunnen ze onafhankelijk van elkaar geflipt worden??

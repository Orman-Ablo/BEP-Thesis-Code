function dataout = cztfuncn(datain,allA,allB,allD,printsize)
% This function evaluates the FT via the n-dimensional czt-algorithm
% arguments: datain = input data, dimensions N1 x N2 x ... x Nn
%            allA,allB,allD = auxiliary vectors computed in prechirpzn, 
%            nx1 cell arrays with vectors of length Nk, Mk, and Lk=Nk+Mk-1
% function value: dataout = output data, dimensions M1 x M2 x ... x Mn
%
% copyright Sjoerd Stallinga, TU Delft, 2023

% retrieve dimensionality, permutation vector for cycling through dimensions
nd = length(size(datain)); 
dimorder = [2:nd,1];

% energy_in = sum(abs(datain(:)).^2); % input energy

% loop over dimensions, apply 1D-CZT to all dimensions of the array
dataout = datain;
for jd = 1:nd
  A = allA{jd};
  B = allB{jd};
  D = allD{jd};
  if printsize
    fprintf('Size cztfuncn input loop %i\n',jd)
    disp(size(dataout))
  end
  dataout = cztfunc1(dataout,A,B,D); 
    % energy_out = sum(abs(dataout(:)).^2);% output energy
    % if energy_in > 0
    % dataout = dataout * sqrt(energy_in / energy_out); % energy conservation
    % end
  if printsize
    fprintf('Size cztfuncn output loop %i\n',jd)
    disp(size(dataout))
  end
  dataout = permute(dataout,dimorder);
  if printsize
    fprintf('Size cztfuncn after permutation loop %i\n',jd)
    disp(size(dataout))
  end
end
  
end
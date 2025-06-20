function [allA,allB,allD] = prechirpzn(allxsize,allqsize,allN,allM)
% This function evaluates the auxiliary vectors for the evaluation of the
% FT via the n-dimensional czt-algorithm
% arguments: xsize = window in real space abs(x)<xsize
%            qsize = window in Fourier space abs(q)<qsize
%            N = # sample points in real space (even)
%            M = # sample points in Fourier space (odd)
% function value: A,B,D = auxiliary vectors of lengths N, M, and L=N+M-1
%
% copyright Sjoerd Stallinga, TU Delft, 2023

% retrieve dimensionality
nd = length(allM); 

% initialization
allA = cell(nd,1);
allB = cell(nd,1);
allD = cell(nd,1);

% loop over dimensions
for jd = 1:nd
  xsize = allxsize(jd);
  qsize = allqsize(jd);
  N = allN(jd);
  M = allM(jd);

  % computation ABD
  L = N+M-1;
  sigma = 2*pi*xsize*qsize/N/M;
  Gfac = (2*xsize/N)*exp(1i*sigma*(1-N)*(1-M));

  allns = 0:(N-1);
  A = exp(2*1i*sigma*allns.*(allns+1-M)); 
  allms = 0:(M-1);
  B = Gfac*exp(2*1i*sigma*allms.*(allms+1-N));
  allnms = 0:max(N,M)+1;
  Vtmp = exp(2*1i*sigma*allnms.^2);

  D = ones(1,L);
  for i=1:M
    D(i) = conj(Vtmp(i));
  end
  for i=M+1:L
    D(i) = conj(Vtmp(L+2-i));
  end
  D = fft(D);

  % store ABD vectors
  allA{jd} = A;
  allB{jd} = B;
  allD{jd} = D;
end

end


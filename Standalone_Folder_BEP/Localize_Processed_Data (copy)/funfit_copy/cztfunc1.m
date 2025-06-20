function dataout = cztfunc1(datain,A,B,D)
% This function evaluates the FT via the czt-algorithm
% arguments: datain = input data, dimensions N1xN2x...xNn
%            A,B,D = auxiliary vectors computed in prechirpz, must have
%            lengths N=N1, M=M1, and L=L1=N1+M1-1
% function value: dataout = output data, dimensions M1xN2x...xNn
%
% copyright Sjoerd Stallinga, TU Delft, 2023

N = length(A);
M = length(B);
L = length(D);

% bring to NxK shape
insizevec = size(datain);
dimsize = prod(insizevec);
K = dimsize/N;



datain = reshape(datain,[N K]);

% make 1D-CZT in first dimension

% % ... below code contains a bug ...
% Amt = repmat(A',1,K);
% Bmt = repmat(B',1,K);
% Dmt = repmat(D',1,K);
% 
% cztin =  zeros(L,K);
% cztin(1:N,:)= Amt.*datain;
% temp = Dmt.*fft(cztin,[],1);
% cztout = ifft(temp,[],1);
% dataout = Bmt.*cztout(1:M,:);

% ... route with transpose to move the FT to column direction works ...
datain = transpose(datain);
Amt = repmat(A,K,1);
Bmt = repmat(B,K,1);
Dmt = repmat(D,K,1);

cztin =  zeros(K,L);
cztin(:,1:N)= Amt.*datain;
temp = Dmt.*fft(cztin,[],2);
cztout = ifft(temp,[],2);
dataout = Bmt.*cztout(:,1:M);
dataout = transpose(dataout);

% bring to M1 x N2 x ... x Nn shape
outsizevec = insizevec;
outsizevec(1) = M;
dataout = reshape(dataout,outsizevec);

end



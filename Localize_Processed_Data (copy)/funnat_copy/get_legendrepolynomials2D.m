function alllegendres2D = get_legendrepolynomials2D(orders,x,y)
% This function computes products of Legendre polynomials in x and y for
% the given set of order combinations. We use products
%
%   P2D(n,m,x,y) = P(n,x)*P(m,y)
%
% where the input is expected as [n1 m1; n2 m2; ....], and where 
%   n,m = 0,1,2,...
%
% copyright Sjoerd Stallinga, TU Delft, 2023

Nleg = size(orders,1); % number of Legendre polynomial products

% change arrays x and y into vector, x and y are assumed to have the same
% array size
xydims = size(x);
Nxy = prod(xydims);
x = reshape(x,[Nxy 1]);
y = reshape(y,[Nxy 1]);

% Evaluation of the Legendre polynomials using the recurrence relation +2,for
% these polynomials
xorders = orders(:,1);
alllegsx = get_legendrepolynomials(xorders,x);
yorders = orders(:,2);
alllegsy = get_legendrepolynomials(yorders,y);
alllegsx = squeeze(alllegsx);
alllegsy = squeeze(alllegsy);

if Nxy == 1 % make sure that first dimension represents nr of points
    alllegsx = alllegsx';
    alllegsy = alllegsy';
end

% Compute the products of the Legendre polynomials in x and y
alllegendres2D = zeros(Nxy,Nleg);
for jleg = 1:Nleg
  alllegendres2D(:,jleg) = alllegsx(:,jleg).*alllegsy(:,jleg);
end

% Reshape computed array of Legendre polynomial products back to original
% array size of x and y
legdims = [xydims,Nleg];
alllegendres2D = reshape(alllegendres2D,legdims);

end


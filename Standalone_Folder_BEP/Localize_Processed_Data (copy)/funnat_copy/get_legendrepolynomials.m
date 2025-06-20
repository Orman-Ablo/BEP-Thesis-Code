function alllegendres = get_legendrepolynomials(orders,x)
% This function computes the Legendre polynomials for the given orders for
% an array of x-values (must be in the range [-1,1].
%
% copyright Sjoerd Stallinga, TU Delft, 2023

Nleg = length(orders);
ordmax = max(orders);

% change array x into a vector
xdims = size(x);
Nx = prod(xdims);
x = reshape(x,[Nx 1]);

% Evaluation of the Legendre polynomials using the recurrence relation for
% these polynomials
legpols = zeros(Nx,ordmax+1);
legpols(:,1) = ones(size(x));
if ordmax>=1
  legpols(:,2) = x;
  if ordmax>=2
    for jord = 3:ordmax+1
      n = jord-2;
      legpols(:,jord) = ((2*n+1)*x.*legpols(:,jord-1)-n*legpols(:,jord-2))/(n+1);
    end
  end
end

% Retain the Legendre polynomials of choice
alllegendres = zeros(Nx,Nleg);
for jleg = 1:Nleg
  jord = orders(jleg)+1;
  alllegendres(:,jleg) = legpols(:,jord);
end

% reshape back to original dimensions of array x
legdims = [xdims,Nleg];
alllegendres = reshape(alllegendres,legdims);


end


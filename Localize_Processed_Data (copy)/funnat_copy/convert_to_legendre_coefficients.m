function [gammas_legendre] = convert_to_legendre_coefficients(gammas_standard)
%convert_to_legendre_coefficients convert standard gamma coefficients to
%legendre coefficients, in order to compare if both lead to the same (CRLB)
%surfaces.

gammas_legendre = zeros(size(gammas_standard));

gammas_legendre(1) = gammas_standard(1) + gammas_standard(4)*(2/3);
gammas_legendre(2) = gammas_standard(2)/sqrt(3);
gammas_legendre(3) = gammas_standard(3)/sqrt(3);
gammas_legendre(4) = gammas_standard(4)*2/(3*sqrt(5));
gammas_legendre(5) = gammas_standard(6);
gammas_legendre(6) = gammas_standard(5);
gammas_legendre(7) = gammas_standard(8)/sqrt(3);
gammas_legendre(8) = gammas_standard(7)/sqrt(3);
gammas_legendre(9) = gammas_standard(9)*(2/3);
gammas_legendre(10) = gammas_standard(11);
gammas_legendre(11) = gammas_standard(10);
gammas_legendre(12) = gammas_standard(12)/sqrt(3);
gammas_legendre(13) = gammas_standard(13);

end
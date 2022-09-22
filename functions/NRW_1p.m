function [eps mu] = NRW_1p(f, GA, GB, ai, t, A, B, epsm, mum)
% f - frequency, Hz;
% S11, S21 - s-parameters measured at the sample faces
% ai = 1/(4*a^2) for a rectangular waveguide, ai = 0 for a coaxial line
% L - thickness of the sample, m
% epsm, mum - permittivity and permeability of the waveguide media

c = 299792458;
k3 = 2*pi*sqrt(epsm.*mum./((c./f).^2)-ai);
Z3 = 2*pi*f.*mum./k3./c;
ZA = Z3.*(1+GA)./(1-GA);
ZB = Z3.*(1+GB)./(1-GB);
c1 = (ZA.*A.*(B-ZB)+ZB.*B.*(ZA-A))./(B-A+ZA-ZB);
U1 = sqrt(c1).*(A-ZA)./(1i*(ZA.*A-c1)); % !!!
%%%
%U2 = -sqrt(c1).*(A-ZA)./(1i*(ZA.*A-c1));
%U = abs(real(U1))+1i*imag(U1);
%Iu = real(U1)>0;
%U = U2;
%U(Iu) = U1(Iu);
%clear Ia U1 U2
%%%

UT = (-1-1i*U1)./(1i.*U1-1);
V = (log(abs(UT))+1i*unwrap(angle(UT)))./(2i*t);

c2 = (V.^2+4*pi^2*ai).*c^2./(2*pi*f).^2;
mu = sqrt(c1.*c2-c1.*(4*ai*c.^2)./(4*f.^2));
eps = c2./mu;

end
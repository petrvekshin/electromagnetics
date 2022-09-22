function [eps mu] = NRW(f, S11, S21, ai, L, epsm, mum, alpha, n)
% f - frequency [Hz];
% S11, S21 - S-parameters measured at the sample faces
% ai = 1/(4*a^2) for TE10 mode, ai = 0 for TEM mode [m^-2]
% L - the thickness of the sample [m]
% epsm, mum - permittivity and permeability of the transmission medium
% alpha - attenuation constant (attenuation due to metal conductivity)
% n - integer number of wavelengths in the sample at the lowest frequency (n = 0, 1, ...)

l2 = (299792458./f).^2;

N = length(f);

X = 0.5*(S11.^2-S21.^2+1)./S11;

Ga = X+sqrt(X.^2-1);
Gb = X-sqrt(X.^2-1);
Ia = abs(Ga)<1;
G1 = Gb;
G1(Ia) = Ga(Ia);
clear Ia Ga Gb

T = (S11+S21-G1)./(1-(S11+S21).*G1);

angT = angle(1./T);
angT = (angT<=0)*2*pi+angT;
ln1T = log(abs(1./T))+1i*unwrap(angT+n*2*pi);
U = -((ln1T./(2*pi*L)).^2)+ln1T.*alpha/(2*L*pi^2)-alpha.^2/(2*pi)^2;
mu = mum.*sqrt(U).*(1+G1)./((1-G1).*sqrt(epsm.*mum./l2-ai));
%mu = ones(N,1);
eps = l2.*(ai+U)./mu;

end
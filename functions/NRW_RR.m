function [eps mu] = NRW_RR(f, A, B, V, U, ai, L, epsm, mum, alpha, n)
% f - frequency [Hz];
% A, B - Sii measured with the MUT
% V, U - Sii measured without the MUT
% ai = 1/(4*a^2) for TE10 mode, ai = 0 for TEM mode [m^-2]
% L - the thickness of the sample [m]
% epsm, mum - permittivity and permeability of the transmission medium
% alpha - attenuation constant (attenuation due to metal conductivity)
% n - integer number of wavelengths in the sample at the lowest frequency (n = 0, 1, ...)

l2 = (299792458./f).^2;

N = length(f);

Q = A - B + U - V + A.*B.*(U-V) + (A-B).*U.*V;
W = sqrt(-4*(B.*V - A.*U).^2 + Q.^2);

Z = 2*(B.*V - A.*U);

Ga = (-Q+W)./Z;
Gb = (-Q-W)./Z;
Ia = abs(Ga)<1;
G1 = Gb;
G1(Ia) = Ga(Ia);
clear Ia Ga Gb Q W Z

Pa = sqrt((G1-B).*(U.*G1-1)./(B.*G1.*(G1-U)-G1+U));
%Pa = sqrt((G1-A).*(V.*G1-1)./(A.*G1.*(G1-V)-G1+V));
Pb = -Pa;
P1 = zeros(N,1);
if (angle(Pa(1)) < 0)
    P1(1) = Pa(1);
else
    P1(1) = Pb(1);
end
for I = 2:N
    if (abs(angle(Pa(I)/P1(I-1))) < pi/2)
        P1(I) = Pa(I);
    else
        P1(I) = Pb(I);
    end
end

T = P1;
clear P1

angT = angle(1./T);
angT = (angT<=0)*2*pi+angT;
ln1T = log(abs(1./T))+1i*unwrap(angT+n*2*pi);
U = -((ln1T./(2*pi*L)).^2)+ln1T.*alpha/(2*L*pi^2)-alpha.^2/(2*pi)^2;
mu = mum.*sqrt(U).*(1+G1)./((1-G1).*sqrt(epsm.*mum./l2-ai));
%mu = ones(N,1);
eps = l2.*(ai+U)./mu;

end
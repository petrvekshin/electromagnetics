function [V,W] = Reflection(f, d, eps, mu, Zi, Zt)
% Zi and Zt are normalized impedances
c = 299792458;
N = length(d);
Zin = Zt;
W = 2*Zin./(Zin+sqrt(mu(:,1)./eps(:,1)));
for I = 1:N
    k = 2*pi*f.*sqrt(eps(:,I).*mu(:,I))/c;
    Z1 = Zin;
    Z2 = sqrt(mu(:,I)./(eps(:,I)));
    Zin = Z2.*(Z1+1i.*Z2.*tan(k*d(I)))./(Z2+1i.*Z1.*tan(k*d(I)));
    if I < N
        W = W.*(Zin+Z2).*exp(-1i*k*d(I))./(Zin+sqrt(mu(:,I+1)./(eps(:,I+1))));
    end
end

V = (Zin-Zi)./(Zin+Zi);
W = W.*(Zin+Z2).*exp(-1i*k*d(N))./(Zin+Zi);
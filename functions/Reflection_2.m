function [V_n,W_n,V_p,W_p] = Reflection_2(f, d, eps, mu, thi, Zi, ki, Zt, kt)
% Zi and Zt are normalized impedances
c = 299792458;
N = length(d);

tht = asind(ki.*sind(thi)./kt);

Zin_n = Zt./cosd(tht); % normal to the plane of incidence (s-pol.)
Zin_p = Zt.*cosd(tht); % along the plane of incidence (p-pol.)

W_n = 2*Zin_n./(Zin_n+sqrt(mu(:,1)./eps(:,1))...
./cosd(asind(ki.*sind(thi)./(2*pi*f.*sqrt(eps(:,1).*mu(:,1))/c))));

W_p = 2*Zin_p./(Zin_p+sqrt(mu(:,1)./eps(:,1))...
.*cosd(asind(ki.*sind(thi)./(2*pi*f.*sqrt(eps(:,1).*mu(:,1))/c))));


for I = 1:N
    k0 = 2*pi*f.*sqrt(eps(:,I).*mu(:,I))./c;
    costh = cosd(asind(ki.*sind(thi)./k0));
    k = k0.*costh;
    Z1_n = Zin_n;
    Z1_p = Zin_p;
    Z2_n = sqrt(mu(:,I)./(eps(:,I)))./costh;
    Z2_p = sqrt(mu(:,I)./(eps(:,I))).*costh;
    Zin_n = Z2_n.*(Z1_n+1i.*Z2_n.*tan(k*d(I)))./(Z2_n+1i.*Z1_n.*tan(k*d(I)));
    Zin_p = Z2_p.*(Z1_p+1i.*Z2_p.*tan(k*d(I)))./(Z2_p+1i.*Z1_p.*tan(k*d(I)));
    if I < N
        costhf = cosd(asind(ki.*sind(thi)./(2*pi*f.*sqrt(eps(:,I+1).*mu(:,I+1))./c)));
        W_n = W_n.*(Zin_n+Z2_n).*exp(-1i*k*d(I))...
        ./(Zin_n+sqrt(mu(:,I+1)./(eps(:,I+1)))./costhf);
        W_p = W_p.*(Zin_p+Z2_p).*exp(-1i*k*d(I))...
        ./(Zin_p+sqrt(mu(:,I+1)./(eps(:,I+1))).*costhf);
    end
end

V_n = (Zin_n-Zi./cosd(thi))./(Zin_n+Zi./cosd(thi));
V_p = (Zin_p-Zi.*cosd(thi))./(Zin_p+Zi.*cosd(thi));

W_n = W_n.*(Zin_n+Z2_n).*exp(-1i*k*d(N))./(Zin_n+Zi./cosd(thi));
W_p = W_p.*(Zin_p+Z2_p).*exp(-1i*k*d(N))./(Zin_p+Zi.*cosd(thi));
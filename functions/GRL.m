function [S] = GRL(f,E,M,SMUT,A1,A2,d,z1,z2)
% Improved Free-Space S-Parameter Calibration
% by Philip G. Bartley, Jr. and Shelley B. Begley

% f is frequency
% E = [E11 E21 E12 E22] is the empty fixture
% M = [M11 M22] is the fixture with a metal plate
% SMUT = [SMUT11 SMUT21 SMUT12 SMUT22] is the fixture with MUT
% A1, A2 are the antennas
% d is the thickness of the metal plate
% z1, z2 are coordinates of MUT

c0 = 299792458;
mu0 = pi*4e-7;

n = length(f);
epsm = 1.000649*ones(n,1);
mum = ones(n,1);
A21 = exp(-2i*pi*f.*sqrt(epsm.*mum)*d./c0);
A12 = A21;
Gair1 = E(:,1)-A1;
Gair2 = E(:,4)-A2;
Gmet1 = M(:,1)-A1;
Gmet2 = M(:,2)-A2;

a = Gmet2.*(A21.*A12).^2.*(Gmet1-Gair1);
b = Gmet1.*Gmet2.*(A21.*A12).^2-Gair1.*Gmet2.*A21.*A12+Gair2.*Gmet1.*A21.*A12;
c = Gair2.*Gmet1.*A21.*A12;

Ga =(-b+sqrt(b.^2-4.*a.*c))./(2*a);
Gb = (-b-sqrt(b.^2-4.*a.*c))./(2*a);
Ia = abs(Ga)<abs(Gb);
O22 = Gb;
O22(Ia) = Ga(Ia);
clear Ia Ga Gb
T22 = -Gair1./(O22.*(Gmet1.*A21.*A12-Gair1.*A21.*A12)+Gmet1.*A21.*A12);
O21O12 = -Gmet1.*(1+O22);
T21T12 = -Gmet2.*(1+T22);

%(A1.*SMUT(:,4).*T22-A1.*A2.*T22+A1.*T21T12-...
%SMUT(:,1).*SMUT(:,4).*T22+SMUT(:,3).*SMUT(:,2).*T22+SMUT(:,1).*A2.*T22-SMUT(:,1).*T21T12)

S11 = (T22.*(A1.*(SMUT(:,4)-A2)+SMUT(:,1).*(A2-SMUT(:,4))+SMUT(:,2).*SMUT(:,3))+T21T12.*(A1-SMUT(:,1)))./...
(A1.*O22.*SMUT(:,4).*T22-O21O12.*SMUT(:,4).*T22-A1.*O22.*A2.*T22+...
A1.*O22.*T21T12+O21O12.*A2.*T22-O21O12.*T21T12-...
O22.*SMUT(:,1).*SMUT(:,4).*T22+O22.*SMUT(:,3).*SMUT(:,2).*T22+O22.*SMUT(:,1).*A2.*T22-O22.*SMUT(:,1).*T21T12);

S22 = (A1.*O22.*SMUT(:,4)-O21O12.*SMUT(:,4)-A1.*O22.*A2+O21O12.*A2-...
O22.*SMUT(:,1).*SMUT(:,4)+O22.*SMUT(:,3).*SMUT(:,2)+O22.*SMUT(:,1).*A2)./(A1.*O22.*SMUT(:,4).*T22-...
O21O12.*SMUT(:,4).*T22-A1.*O22.*A2.*T22+A1.*O22.*T21T12+O21O12.*A2.*T22-...
O21O12.*T21T12-O22.*SMUT(:,1).*SMUT(:,4).*T22+O22.*SMUT(:,3).*SMUT(:,2).*T22+O22.*SMUT(:,1).*A2.*T22-O22.*SMUT(:,1).*T21T12);

S11 = S11.*exp(4i*pi*f.*sqrt(epsm.*mum)*z1./c0);
S22 = S22.*exp(4i*pi*f.*sqrt(epsm.*mum)*(d-z2)./c0);

O21T12 = E(:,2).*(1-A21.*A12.*O22.*T22)./A21;
O12T21 = E(:,3).*(1-A21.*A12.*O22.*T22)./A12;
DNS21S12 = A1.*O22.*SMUT(:,4).*T22-O21O12.*SMUT(:,4).*T22-...
A1.*O22.*A2.*T22+A1.*O22.*T21T12+O21O12.*A2.*T22-O21O12.*T21T12-...
O22.*SMUT(:,1).*SMUT(:,4).*T22+O22.*SMUT(:,3).*SMUT(:,2).*T22+O22.*SMUT(:,1).*A2.*T22-O22.*SMUT(:,1).*T21T12;

S21 = -O12T21.*SMUT(:,2)./DNS21S12;
S12 = -O21T12.*SMUT(:,3)./DNS21S12;

S21 = S21.*exp(-2i*pi*f.*sqrt(epsm.*mum)*(z2-z1-d)./c0);
S12 = S12.*exp(-2i*pi*f.*sqrt(epsm.*mum)*(z2-z1-d)./c0);

S = [S11,S21,S12,S22];
end


function [mRCS] = fast_mRCS_1x_convex_with_a_covering(M,k0,H,lambda,d,eps,mu,Zt,kt)
% Fast mRCS computation of a convex PEC body using PO with numerical integration

% Wavenumber
wn = 2*pi./lambda;

% Finding visible polygons
V1 = M(:,13).*k0(1,1)+M(:,14).*k0(1,2)+M(:,15).*k0(1,3);
V = V1<-0.0001;
%m = sum(V);
V1 = V1(V,:);
M = M(V,:);
thi = acosd(-V1);

Hpp = incwave2twovectors(k0, H, M(:,13:15));
L1 = sqrt(Hpp(:,1).^2+Hpp(:,2).^2+Hpp(:,3).^2);
%L2 = sqrt(Hpp(:,4).^2+Hpp(:,5).^2+Hpp(:,6).^2);
clear Hpp
alpha = acosd(L1);
clear L1

[V_n,V_p] = Reflection_3(lambda, thi, d, eps, mu, Zt, kt);
H1 = V_n.*cosd(alpha).^2;
H2 = V_p.*sind(alpha).^2;
clear V_n V_p
HS = H1+H2;
clear H1 H2

Hf = -M(:,16).*(k0(1,1).*M(:,13)...
  +k0(1,2).*M(:,14)+k0(1,3).*M(:,15)).*exp(-2i.*wn.*(k0(1,1).*M(:,10)+...
  k0(1,2).*M(:,11)+k0(1,3).*M(:,12))).*HS;

% H-field computation
Hf = sum(Hf);

% mRCS computation
mRCS = (4*pi/lambda^2)*abs(Hf)^2;

end
% Measurements of planar microwave circuits
% using an improved TRL calibration method
% Y. Liu, L. Tong, Y. Tian, and B. Gao
% Progress In Electromagnetics Research, Vol. 109, 263�278, 2010

clear all
close all
clc

% Line
[SL,freq] = S_Parameters_Loading('../TRL_Line.s2p');
n = length(freq);

% Thru
ST = S_Parameters_Loading('../TRL_Thru.s2p');

% Reflect 1
SR1 = S_Parameters_Loading('../TRL_Short_1.s2p');

% Reflect 2
SR2 = S_Parameters_Loading('../TRL_Short_2.s2p');

% Measurement
SM = S_Parameters_Loading('../Material_wo_Calibration.s2p');

b = ST(:,2).*ST(:,3)+SL(:,3).*SL(:,2)+(SL(:,4)-ST(:,4)).*(ST(:,1)-SL(:,1));

p1X = -b./(2*(-ST(:,2)).*SL(:,3));
p2X = sqrt(b.^2-4*ST(:,2).*SL(:,3).*SL(:,2).*ST(:,3))./(2*(-ST(:,2)).*SL(:,3));

J = angle(p1X+p2X) < 0;
X = zeros(n,1);
X(J) = p1X(J)+p2X(J);
X(~J) = p1X(~J)-p2X(~J);

clear p1X p2X J

A = (ST(:,4)-SL(:,4))./(ST(:,3)-SL(:,2).*X);

W = ((ST(:,4)-SL(:,4))./(ST(:,3)-SL(:,3).*X)).*...
((SR1(:,1)-ST(:,1))./(ST(:,2).*(1-A)))+(A./(1-A));

V = ((ST(:,1)-SL(:,1))./(ST(:,2)-SL(:,2).*X)).*...
((SR2(:,4)-ST(:,4))./(ST(:,3).*(1-A)))+(A./(1-A));

G1 = sqrt(W.*V./((1+W).*(1+V).*A));
J = angle(G1) < 0;
G = zeros(n,1);
G(J) = G1(J);
G(~J) = -G1(~J);
clear G1 J

S22A = W./(G.*(1+W));
S11B = V./(G.*(1+V));
S11A = ST(:,1)-((1-A.*(X.^2)).*(ST(:,1)-SL(:,1))./(1-(X.^2)));
S22B = ST(:,4)-((1-A.*(X.^2)).*(ST(:,4)-SL(:,4))./(1-(X.^2)));
T = ST(:,2).*(1-A);
P = ST(:,3).*(1-A);
Z = (ST(:,1)-SL(:,1))./(S11B.*((1./(1-A))-((X.^2)./(1-A.*(X.^2)))));
Y = (ST(:,4)-SL(:,4))./(S22A.*((1./(1-A))-((X.^2)./(1-A.*(X.^2)))));
B = (1+((SM(:,1)-S11A)./Z).*S22A).*(1+((SM(:,4)-S22B)./Y).*S11B)-...
(SM(:,2)./T).*(SM(:,3)./P).*S11B.*S22A;

S21X = SM(:,2)./(T.*B);
S12X = SM(:,3)./(P.*B);

plot(freq,angle(S12X))






























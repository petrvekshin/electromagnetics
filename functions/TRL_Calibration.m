clear all
close all
clc
c0 = 299792458;

eps = 1.000649;
fc1 = 7e9;
fc2 = 13e9;
f1 = 8.5e9;
f2 = 12.4e9;
awg = [22:0.01:24]'*1e-3;
N = length(awg);

% Line
[MD,freq] = S_Parameters_Loading('../L.s2p');
fcI = freq >= fc1;
fcI = logical(fcI.*(freq <= fc2));
ncI = length(fcI);
freq = freq(fcI);
n = length(freq);
fI = freq >= f1;
fI = logical(fI.*(freq <= f2));
nI = length(fI);

RD = Spar2Rpar(MD(fcI,:));
clear MD

%W = exp(-1i*sqrt(eps*(2*pi*freq/c).^2-(pi/awg)^2)*0.00659);

% Thru
RT = Spar2Rpar(S_Parameters_Loading('../T.s2p'));
RT = RT(fcI,:);

% Reflect 1, 2
MS1 = S_Parameters_Loading('../R.s2p');
w1 = MS1(fcI,1);
w2 = MS1(fcI,4);
clear MS1

% Reflect 2
%MS2 = S_Parameters_Loading('../Reflect1.s2p');
%w2 = MS2(fcI,4);
%clear MS2

% Measurement
%RM = Spar2Rpar(S_Parameters_Loading('../Line.s2p'));
%RM = RM(fcI,:);
RM = RD;

invRT = invM(RT);
T = [RD(:,1).*invRT(:,1)+RD(:,3).*invRT(:,2),...
    RD(:,2).*invRT(:,1)+RD(:,4).*invRT(:,2),...
    RD(:,1).*invRT(:,3)+RD(:,3).*invRT(:,4),...
    RD(:,2).*invRT(:,3)+RD(:,4).*invRT(:,4)];

D = (T(:,4)-T(:,1)).^2+4*T(:,2).*T(:,3);
x1 = (-(T(:,4)-T(:,1))+sqrt(D))./(2*T(:,2));
x2 = (-(T(:,4)-T(:,1))-sqrt(D))./(2*T(:,2));
clear D

I = abs(x1)>abs(x2);

a_c = zeros(n,1);
a_c(I) = x1(I);
a_c(~I) = x2(~I);

b = zeros(n,1);
b(I) = x2(I);
b(~I) = x1(~I);

clear x1 x2 I

e2gL = (T(:,2).*b+T(:,4))./(T(:,2).*a_c+T(:,4));

d = RT(:,1)./RT(:,4);
f = RT(:,2)./RT(:,4);
e = RT(:,3)./RT(:,4);

r22rho22 = RT(:,4).*(1-e./a_c)./(1-b./a_c);

g = (f-d./a_c)./(1-e./a_c);

beta_alpha = (e-b)./(d-b.*f);

alpha_a = (d-b.*f)./(1-e./a_c);

a = sqrt((w1-b).*(1+w2.*beta_alpha).*(d-b.*f)./...
((w2+g).*(1-w1./a_c).*(1-e./a_c)));
c1 = a./a_c;
G = (w1-b)./(a-c1.*w1);

% Assuming that G is close to -1
J = real(G)<0;
a(~J) = -a(~J);
G(~J) = -G(~J);

alpha = alpha_a./a;

% Using the known reflection coefficient that is -1
%a0 = (w1-b)./(w1./a_c-1);
%alpha0 = (w2+g)./(-w2.*beta_alpha-1);

c = a./a_c;

beta = beta_alpha.*alpha;

D_DUT = (1./r22rho22).*(1./alpha_a).*(1./(1-b./a_c)).*(1./(1-g.*beta_alpha));

R11 = RM(:,1)-g.*(RM(:,3)-RM(:,4).*b)-RM(:,2).*b;
R21 = RM(:,2).*a-RM(:,1).*c-g.*(RM(:,4).*a-RM(:,3).*c);
R12 = alpha.*(RM(:,3)-RM(:,4).*b)-beta.*(RM(:,1)-RM(:,2).*b);
R22 = alpha.*(RM(:,4).*a-RM(:,3).*c)-beta.*(RM(:,2).*a-RM(:,1).*c);

R = [R11.*D_DUT R21.*D_DUT R12.*D_DUT R22.*D_DUT];

S = Rpar2Spar(R);

FontSize = 14;

%f1 = figure('Name','Scattering parameters');
%subplot(2,2,1)
%plot(freq/10^9,20*log10(abs(S(:,1))),'r','LineWidth',1.5);
%hold on
%plot(freq/10^9,20*log10(abs(S(:,4))),'b','LineWidth',1.5);
%title('\rm\itS\rm_{11}, \itS\rm_{22}','FontSize',14);
%xlabel('Frequency, GHz','FontSize',FontSize);
%ylabel('Magnitude, dB','FontSize',FontSize);
%grid on
%set(gca,'FontSize',FontSize)

%subplot(2,2,2)
%plot(freq/10^9,angle(S(:,1))*180/pi,'r','LineWidth',1.5);
%hold on
%plot(freq/10^9,angle(S(:,4))*180/pi,'b','LineWidth',1.5);
%title('\rm\itS\rm_{11}, \itS\rm_{22}','FontSize',14);
%xlabel('Frequency, GHz','FontSize',FontSize);
%ylabel('Phase, deg','FontSize',FontSize);
%grid on
%set(gca,'FontSize',FontSize)

%subplot(2,2,3)
%plot(freq/10^9,20*log10(abs(S(:,2))),'r','LineWidth',1.5);
%hold on
%plot(freq/10^9,20*log10(abs(S(:,3))),'b','LineWidth',1.5);
%title('\rm\itS\rm_{21}, \itS\rm_{12}','FontSize',14);
%xlabel('Frequency, GHz','FontSize',FontSize);
%ylabel('Magnitude, dB','FontSize',FontSize);
%grid on
%set(gca,'FontSize',FontSize)

%subplot(2,2,4)
%plot(freq/10^9,unwrap(angle(S(:,2)))*180/pi,'r','LineWidth',1.5);
%hold on
%plot(freq/10^9,unwrap(angle(S(:,3)))*180/pi,'b','LineWidth',1.5);
%hold on
%plot(freq/10^9,unwrap(angle(W))*180/pi,'g','LineWidth',2);
%title('\rm\itS\rm_{21}, \itS\rm_{12}','FontSize',14);
%xlabel('Frequency, GHz','FontSize',FontSize);
%ylabel('Phase, deg','FontSize',FontSize);
%grid on
%set(gca,'FontSize',FontSize)

S11 = S(:,1);
S21 = S(:,2);
S12 = S(:,3);
S22 = S(:,4);
f = freq;

A21 = unwrap(-angle(S21));
A12 = unwrap(-angle(S12));
RS21 = zeros(N,1);
RM21 = zeros(N,1);
RS12 = zeros(N,1);
RM12 = zeros(N,1);
d21 = zeros(N,1);
d12 = zeros(N,1);

for I = 1:N
    l_g = 2*pi./sqrt(eps*(2*pi*freq/c0).^2-(pi/awg(I)).^2);
    D21 = A21.*l_g/2/pi;
    D12 = A12.*l_g/2/pi;
    d21(I) = median(D21(fI));
    d12(I) = median(D12(fI));
    phD21 = 2*pi*d21(I)./l_g;
    phD12 = 2*pi*d12(I)./l_g;
    RS21(I,1) = sum(abs(phD21(fI)-A21(fI)))/n;
    RS12(I,1) = sum(abs(phD12(fI)-A12(fI)))/n;
    RM21(I,1) = median(abs(phD21(fI)-A21(fI)));
    RM12(I,1) = median(abs(phD12(fI)-A12(fI)));
end
%[i1,i2] = min(R);
%[i3,i4] = min(R2);

% Plot properties
FontSize = 22;
AxisColor = [.95 .95 .95];

figure('Position',[11 51 1800 900],'Name','Error','Color',[0 0 0]);

subplot(2,1,1)
plot(awg*1000,d21*1000,'LineWidth',2.5,'Color',[1 .3 .3]);
hold on
plot(awg*1000,d12*1000,'LineWidth',2.5,'Color',[.3 .3 1]);
grid on
xlabel('a, mm','FontSize',FontSize,'Color',AxisColor);
ylabel('Length, mm','FontSize',FontSize,'Color',AxisColor);
xlim([awg(1),awg(end)]*1000)
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,1,2)
plot(awg*1000,RS21*180/pi,'LineWidth',2.5,'Color',[.3 .3 1]);
hold on
plot(awg*1000,RM21*180/pi,'LineWidth',2.5,'Color',[1 .3 .3]);
hold on
plot(awg*1000,RS12*180/pi,'LineWidth',2.5,'Color',[.2 .2 0.7]);
hold on
plot(awg*1000,RM12*180/pi,'LineWidth',2.5,'Color',[0.7 .2 .2]);
grid on
xlabel('a, mm','FontSize',FontSize,'Color',AxisColor);
ylabel('Error, deg','FontSize',FontSize,'Color',AxisColor);
xlim([awg(1),awg(end)]*1000)

set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)

%legend({strcat('Mean',' (',num2str(awg(i2)*1000,4),',',num2str(d(i2)*1000,3),')'),...
%strcat('Median',' (',num2str(awg(i4)*1000,4),',',num2str(d(i4)*1000,3),')')},...
%'location', 'northwest')

wcck = [23.015];
disp([wcck,(interp1(awg*1000,d21*1000,wcck)+interp1(awg*1000,d12*1000,wcck))/2])
IERM = zeros(4,2);
[IERM(1,1),IERM(1,2)] = min(RM21);
[IERM(2,1),IERM(2,2)] = min(RM12);
[IERM(3,1),IERM(3,2)] = min(RS21);
[IERM(4,1),IERM(4,2)] = min(RS12);
disp(awg(IERM(:,2))*1000)

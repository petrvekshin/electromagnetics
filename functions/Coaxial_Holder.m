clear all
close all
clc

c = 299792458;
epsa = 1.000649;
f1 = 1e9;
f2 = 18e9;

[S,f] = S_Parameters_Loading('.../Coax_Holder3.s2p');
fI = f >= f1;
fI = logical(fI.*(f <= f2));
f = f(fI);
n = length(f);

S11 = S(fI,1);
S21 = S(fI,2);
S12 = S(fI,3);
S22 = S(fI,4);
clear S

A = unwrap(-angle(S21));
A2 = unwrap(-angle(S12));
l_g = 2*pi./sqrt(epsa*(2*pi*f/c).^2);
D = A.*l_g/2/pi;
D2 = A2.*l_g/2/pi;
d = median(D);
d2 = median(D2);
disp([d*1000 d*sqrt(epsa)/c*1e12])
disp([d2*1000 d*sqrt(epsa)/c*1e12])

f1 = figure('Name','S21');
plot(f,A*180/pi,'LineWidth',2,'r');
hold on
plot(f,(2*pi*d./l_g)*180/pi,'LineWidth',2,'b');
grid on
xlabel('Frequency, GHz','FontSize',14);
ylabel('Phase, deg','FontSize',14);

f2 = figure('Name','S12');
plot(f,A2*180/pi,'LineWidth',2,'r');
hold on
plot(f,(2*pi*d2./l_g)*180/pi,'LineWidth',2,'b');
grid on
xlabel('Frequency, GHz','FontSize',14);
ylabel('Phase, deg','FontSize',14);

OL1 = median(-2*log(abs(S21))*(50/(107.21*1e-12))./sqrt(f/1e9));

OL = 1.0e9;
ILdB = OL*10*107.21*1e-12.*sqrt(f./1e9)./(log(10)*50);
f3 = figure('Name','S21, S12');
plot(f,20*log10(abs(S21)),'LineWidth',2,'r');
hold on
plot(f,20*log10(abs(S12)),'LineWidth',2,'b');
hold on
plot(f,-ILdB,'LineWidth',2,'k');
grid on
xlabel('Frequency, GHz','FontSize',14);
ylabel('Magnitude, dB','FontSize',14);

clear all
close all
clc

c = 299792458;
epsa = 1.000649;
a = [22.6:0.01:23.4]'*1e-3;
N = length(a);
%P = 0.8;
f1 = 8.0e9;
f2 = 12.8e9;

L0 =  6.615e-3;
dL0 =  .020e-3;
A0 = 23.018e-3;
dA0 =  .070e-3;
L = [L0-dL0, L0+dL0];
A = [A0-dA0, A0+dA0]; 

[S,f] = S_Parameters_Loading('../sample.s2p');
fI = f >= f1;
fI = logical(fI.*(f <= f2));
f = f(fI);
n = length(f);

S11 = S(fI,1);
S21 = S(fI,2);
S12 = S(fI,3);
S22 = S(fI,4);
clear S

A21 = unwrap(-angle(S21));
A12 = unwrap(-angle(S12));
RS21 = zeros(N,1);
RM21 = zeros(N,1);
RS12 = zeros(N,1);
RM12 = zeros(N,1);
d21 = zeros(N,1);
d12 = zeros(N,1);

for I = 1:N
    l_g = 2*pi./sqrt(epsa*(2*pi*f/c).^2-(pi/a(I)).^2);
    D21 = A21.*l_g/2/pi;
    D12 = A12.*l_g/2/pi;
    d21(I) = median(D21);
    d12(I) = median(D12);
    phD21 = 2*pi*d21(I)./l_g;
    phD12 = 2*pi*d12(I)./l_g;
    RS21(I,1) = sum(abs(phD21-A21))/n;
    RS12(I,1) = sum(abs(phD12-A12))/n;
    RM21(I,1) = median(abs(phD21-A21));
    RM12(I,1) = median(abs(phD12-A12));
end

L1 = interp1(a,d21,23.018e-3);
L2 = interp1(a,d12,23.018e-3);

if [L1 >= L(1)]*[L1 <= L(2)]*[L2 >= L(1)]*[L2 <= L(2)] == 1
  String1 = 'Length is OK';
else
  String1 = 'Length is NOT OK';
end

[is12,js12] = min(RS12);
[is21,js21] = min(RS21);
[im12,jm12] = min(RM12);
[im21,jm21] = min(RM21);

if [a(min([js12,js21,jm12,jm21])) >= A(1)]*[a(max([js12,js21,jm12,jm21])) <= A(2)] == 1
  String2 = 'Wide side is OK';
else
  String2 = 'Wide side is NOT OK';
end

disp(String1);
disp(String2);

% Plot properties
FontSize = 22;
AxisColor = [.95 .95 .95];

f1 = figure('Position',[11 51 1660 884],'Name','Error','Color',[0 0 0]);

subplot(2,1,1)
plot(a*1000,d12*1000,'LineWidth',2.5,'Color',[.3 .3 1]);
hold on
plot(a*1000,d21*1000,'LineWidth',2.5,'Color',[1 .3 .3]);
title(String1,'FontSize',FontSize,'Color',AxisColor);
grid on
xlabel('a, mm','FontSize',FontSize,'Color',AxisColor);
ylabel('Length, mm','FontSize',FontSize,'Color',AxisColor);
xlim([a(1),a(end)]*1000)
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,1,2)
plot(a*1000,RS12*180/pi,'LineWidth',2.5,'Color',[.3 .3 1]);
hold on
plot(a*1000,RM12*180/pi,'LineWidth',2.5,'Color',[1 .3 .3]);
hold on
plot(a*1000,RS21*180/pi,'LineWidth',2.5,'Color',[.2 .2 0.7]);
hold on
plot(a*1000,RM21*180/pi,'LineWidth',2.5,'Color',[0.7 .2 .2]);
title(String2,'FontSize',FontSize,'Color',AxisColor);
grid on
xlabel('a, mm','FontSize',FontSize,'Color',AxisColor);
ylabel('Error, deg','FontSize',FontSize,'Color',AxisColor);
xlim([a(1),a(end)]*1000)

set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)

%legend({strcat('Mean',' (',num2str(a(i2)*1000,4),',',num2str(d(i2)*1000,3),')'),...
%strcat('Median',' (',num2str(a(i4)*1000,4),',',num2str(d(i4)*1000,3),')')},...
%'location', 'northwest')

%subplot(2,1,2)
%plot(f,A*180/pi,'LineWidth',2,'r');
%hold on
%plot(f,phD*180/pi,'LineWidth',2,'b');
%grid on
%xlabel('Frequency, GHz','FontSize',14);
%ylabel('Phase, deg','FontSize',14);

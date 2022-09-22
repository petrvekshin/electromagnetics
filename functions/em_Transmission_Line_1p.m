clear all
close all
clc

file1 = '../GA_FP.s1p';
file2 = '../GB_FP.s1p';

[filePath1,fileName1,fileExt1] = fileparts(file1);
[filePath2,fileName2,fileExt2] = fileparts(file2);

TL = 'WG';          % Type of transmission line (WG or CA)
t = 0.00294;
d1 = 0;
d2 = 0.00659-t;
L1 = 0.00659-t;
L2 = 0;

% Waveguide parameters
a = .02302;         % The wide side of the waveguide
b = .01004;         % The short side the waveguide
d = .00999;         % The short side of the sample

% Coaxial parameters
Dout  = .00700;     % The outer diameter of the coaxial line
Din   = .00304;     % The inner diameter of the coaxial line
Douts = .00696;     % The outer diameter of the sample
Dins  = .00304;     % The inner diameter of the sample

freq1 = 8e9;
freq2 = 12e9;

c = 299792458;
mu0 = 4*pi*10^-7; 

[filePath1,fileName1,fileExt1] = fileparts(file1);
[filePath2,fileName2,fileExt2] = fileparts(file2);
S1 = importdata(file1,' ');
GA = S1(:,2)+1i*S1(:,3);
S2 = importdata(file2,' ');
GB = S2(:,2)+1i*S2(:,3);
f = S1(:,1);
n = length(f);
clear S1 S2

% The media of transmission line
epsm = 1.000649*ones(n,1);
mum = ones(n,1);

if TL == 'WG'
    ai = 1/(4*a^2);
else
    ai = 0;
end

GA = GA.*exp((2i*pi*sqrt(epsm.*mum./((c./f).^2)-ai))*2*L1);
GB = GB.*exp((2i*pi*sqrt(epsm.*mum./((c./f).^2)-ai))*2*L2);

A = 1i*sqrt(mum./epsm).*tan(2*pi*sqrt(epsm.*mum./((c./f).^2)-ai)*d1);
B = 1i*sqrt(mum./epsm).*tan(2*pi*sqrt(epsm.*mum./((c./f).^2)-ai)*d2);
[eps1,mu1] = NRW_1p(f, GA, GB, ai, t, A, B, epsm, mum, epsm, mum);

% Air gap correction
if TL == 'WG'
    dlt = (b-d)/b;
    eps1_c = eps1*(1-dlt)./(1-eps1*dlt./epsm);
    mu1_c = (mu1-dlt*mum)/(1-dlt);
else
    % This may be incorrect
    L1 = log(Dins/Din)+log(Dout/Douts);
    L2 = log(Douts/Dins);
    L3 = log(Dout/Din);
    eps1_cr = real(eps1)*L2./(L3-real(eps1)*L1);
    eps1_c = eps1_cr+1i*(eps1_cr.*imag(eps1)./real(eps1))*L3./(L3-L1*real(eps1).*(1+(-imag(eps1)./real(eps1)).^2));
    mu1_c = real(mu1)*(L3-L1)/L2+1i*imag(mu1)*L3/L2;
end

% Plot properties
FontSize = 20;
AxisColor = [.95 .95 .95];
LineWidth = 1.5;
Line1Color = [1 .2 .2];
Line2Color = [.2 .2 1];
DarkLine = .6; % Color coefficients of dark lines. DarkLineColor = LineColor*DarkLine

f1 = figure('Position',[11 51 948 900],'Name','Scattering parameters','Color',[0 0 0]);
subplot(1,2,1)
plot(f/1e9,20*log10(abs(GB)),'-','Color',[0 0 1],'LineWidth',2);
hold on
plot(f/1e9,20*log10(abs(GA)),'-','Color',[1 0 0],'LineWidth',2);
title('\rm\itS\rm_{11}, \itS\rm_{22}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
ylabel('Magnitude, dB','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)

subplot(1,2,2)
plot(f/1e9,angle(GB)*180/pi,'-','Color',[0 0 1],'LineWidth',2);
hold on
plot(f/1e9,angle(GA)*180/pi,'-','Color',[1 0 0],'LineWidth',2);
ylim([-180 180]);
title('\rm\itS\rm_{11}, \itS\rm_{22}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
ylabel('Phase, deg','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'YTick',-180:90:180,'FontSize',FontSize,...
'Color',[0 0 0],'YColor',AxisColor,'XColor',AxisColor);

f2 = figure('Position',[963 51 948 900],'Name','Permittivity & permeability','Color',[0 0 0]);
subplot(2,2,1)
plot(f/1e9,real(eps1),'-','Color',Line1Color*DarkLine,'LineWidth',LineWidth);
hold on
plot(f/1e9,real(eps1_c),'-','Color',Line1Color,'LineWidth',LineWidth);
title('\epsilon_r^{\prime}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)
ylim([4 6])

subplot(2,2,2)
plot(f/1e9,-imag(eps1),'-','Color',Line1Color*DarkLine,'LineWidth',LineWidth);
hold on
plot(f/1e9,-imag(eps1_c),'-','Color',Line1Color,'LineWidth',LineWidth);
title('\epsilon_r^{\prime\prime}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,2,3)
plot(f/1e9,real(mu1),'-','Color',Line1Color*DarkLine,'LineWidth',LineWidth);
hold on
plot(f/1e9,real(mu1_c),'-','Color',Line1Color,'LineWidth',LineWidth);
title('\mu_r^{\prime}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,2,4)
plot(f/1e9,-imag(mu1),'-','Color',Line1Color*DarkLine,'LineWidth',LineWidth);
hold on
plot(f/1e9,-imag(mu1_c),'-','Color',Line1Color,'LineWidth',LineWidth);
title('\mu_r^{\prime\prime}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)



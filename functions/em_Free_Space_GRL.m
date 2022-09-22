% Free-Space
% GRL Calibration
clear all
close all
clc

% Fixture with the MUT (*.s2p)
fileM = '../MUT.s2p';
% MUT position [m]
z = [-.2, 3.1]*1e-3;

% Fixture with the metal plate (*.s2p)
fileP = '../Metal_plate.s2p';
% The thickness of the metal plate [m]
d = 3e-3;

% Empty fixture (*.s2p)
fileE = '../Empty_fixture.s2p';

% Antenna 1 with gating (*.s1p)
fileO = '../Antenna1.s1p';
% Antenna 2 with gating (*.s1p)
fileT = '../Antenna2.s1p';

% Integer number of wavelengths in the sample at the lowest frequency (0, 1, ...)
NRW_root = 0;

% Loss data format: 'tan' for loss tangent
lossf = 'im';

% A degree of a polynomial (eps', eps'', mu', mu'')
PFN = [1, 1, 1, 1];

% Plot properties
FontSize = 20;
AxisColor = [160, 160, 160]./255;
LineWidth1 = 1.5;                   % S11, S21, S12, S22
LineWidth2 = 3;                     % fit
MarkerSize = 6;                     % eps, mu
Color0 = [0, 0, 0];                 % background
Color1 = [.7, .2, .2];              % S11, S21
Color2 = [.2, .2, .7];              % S22, S12          
Color3 = [255, 245, 0]./255;        % fit

[MUT,f] = S_Parameters_Loading(fileM);
P = S_Parameters_Loading(fileP);
P(:,2:3) = [];
E = S_Parameters_Loading(fileE);
A1 = S_Parameters_Loading(fileO);
A2 = S_Parameters_Loading(fileT);

S = GRL(f,E,P,MUT,A1,A2,d,z(1),z(2));
S11 = S(:,1);
S21 = S(:,2);
S12 = S(:,3);
S22 = S(:,4);

n = length(f);
epsm = 1.000649*ones(n,1);
mum = ones(n,1);
alpham = zeros(n,1);
[eps1,mu1] = NRW(f,S11,S21,0,z(2)-z(1),epsm,mum,alpham,NRW_root);
[eps2,mu2] = NRW(f,S22,S12,0,z(2)-z(1),epsm,mum,alpham,NRW_root);

PF1 = polyfit(f,real(eps1+eps2)./2,PFN(1));
PF2 = polyfit(f,-imag(eps1+eps2)./2,PFN(2));
PF3 = polyfit(f,real(mu1+mu2)./2,PFN(3));
PF4 = polyfit(f,-imag(mu1+mu2)./2,PFN(4));

eps_fit = polyval(PF1,f)-1i*polyval(PF2,f).*(polyval(PF2,f)>0);
mu_fit = polyval(PF3,f)-1i*polyval(PF4,f).*(polyval(PF4,f)>0);

% Scattering parameters
f1 = figure('Units', 'Normalized', 'Position', [.005 .09 .99 .79],'Name','S-parameters','Color',Color0);
subplot(2,2,1)
plot(f/1e9,20*log10(abs(S22)),'-','Color',Color2,'LineWidth',LineWidth1);
hold on
plot(f/1e9,20*log10(abs(S11)),'-','Color',Color1,'LineWidth',LineWidth1);
title('\rm\itS\rm_{11}, \itS\rm_{22}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
ylabel('Magnitude, dB','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,2,2)
plot(f/1e9,angle(S22)*180/pi,'-','Color',Color2,'LineWidth',LineWidth1);
hold on
plot(f/1e9,angle(S11)*180/pi,'-','Color',Color1,'LineWidth',LineWidth1);
ylim([-180 180]);
title('\rm\itS\rm_{11}, \itS\rm_{22}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
ylabel('Phase, deg','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'YTick',-180:90:180,'FontSize',FontSize,...
'Color',Color0,'YColor',AxisColor,'XColor',AxisColor);

subplot(2,2,3)
plot(f/1e9,20*log10(abs(S12)),'-','Color',Color2,'LineWidth',LineWidth1);
hold on
plot(f/1e9,20*log10(abs(S21)),'-','Color',Color1,'LineWidth',LineWidth1);
title('\rm\itS\rm_{21}, \itS\rm_{12}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
ylabel('Magnitude, dB','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,2,4)
plot(f/1e9,angle(S12)*180/pi,'-','Color',Color2,'LineWidth',LineWidth1);
hold on
plot(f/1e9,angle(S21)*180/pi,'-','Color',Color1,'LineWidth',LineWidth1);
ylim([-180 180]);
title('\rm\itS\rm_{21}, \itS\rm_{12}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
ylabel('Phase, deg','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'YTick',-180:90:180,'FontSize',FontSize,...
'Color',Color0,'YColor',AxisColor,'XColor',AxisColor);

% Permittivity & permeability
f2 = figure('Units', 'Normalized', 'Position', [.005 .09 .99 .79],'Name','Permittivity & permeability','Color',Color0);
subplot(2,2,1)
plot(f/1e9,real(eps2),'.','Color',Color2,'MarkerSize',MarkerSize);                      % eps'_2
hold on
plot(f/1e9,real(eps1),'.','Color',Color1,'MarkerSize',MarkerSize);                      % eps'_1
hold on
plot(f/1e9,real(eps_fit),'-','Color',Color3,'LineWidth',LineWidth2);                    % eps'_fit
title('\epsilon_r^\prime','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,2,2)
if isequal(lossf,'tan') == 1
plot(f/1e9,-imag(eps2)./real(eps2),'.','Color',Color2,'MarkerSize',MarkerSize);         % tan_e_2
hold on
plot(f/1e9,-imag(eps1)./real(eps1),'.','Color',Color1,'MarkerSize',MarkerSize);         % tan_e_1
hold on
plot(f/1e9,-imag(eps_fit)./real(eps_fit),'-','Color',Color3,'LineWidth',LineWidth2);    % tan_e_fit
title('tan\delta_e','FontSize',FontSize,'Color',AxisColor);
else
plot(f/1e9,-imag(eps2),'.','Color',Color2,'MarkerSize',MarkerSize);                     % eps''_2
hold on
plot(f/1e9,-imag(eps1),'.','Color',Color1,'MarkerSize',MarkerSize);                     % eps''_1
hold on
plot(f/1e9,-imag(eps_fit),'-','Color',Color3,'LineWidth',LineWidth2);                   % eps''_fit
title('\epsilon_r^{\prime\prime}','FontSize',FontSize,'Color',AxisColor);
end
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,2,3)
plot(f/1e9,real(mu2),'.','Color',Color2,'MarkerSize',MarkerSize);                       % mu'_2
hold on
plot(f/1e9,real(mu1),'.','Color',Color1,'MarkerSize',MarkerSize);                       % mu'_1
hold on
plot(f/1e9,real(mu_fit),'-','Color',Color3,'LineWidth',LineWidth2);                     % mu'_fit
title('\mu_r^{\prime}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,2,4)
if isequal(lossf,'tan') == 1
plot(f/1e9,-imag(mu2)./real(mu2),'.','Color',Color2,'MarkerSize',MarkerSize);           % tan_m_2
hold on
plot(f/1e9,-imag(mu1)./real(mu1),'.','Color',Color1,'MarkerSize',MarkerSize);           % tan_m_1
hold on
plot(f/1e9,-imag(mu_fit)./real(mu_fit),'-','Color',Color3,'LineWidth',LineWidth2);      % tan_m_fit
title('tan\delta_m','FontSize',FontSize,'Color',AxisColor);
else
plot(f/1e9,-imag(mu2),'.','Color',Color2,'MarkerSize',MarkerSize);                      % mu''_2
hold on
plot(f/1e9,-imag(mu1),'.','Color',Color1,'MarkerSize',MarkerSize);                      % mu''_1
hold on
plot(f/1e9,-imag(mu_fit),'-','Color',Color3,'LineWidth',LineWidth2);                    % mu''_fit
title('\mu_r^{\prime\prime}','FontSize',FontSize,'Color',AxisColor);
end
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)

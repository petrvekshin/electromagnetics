% Air-filled transmission line (rectangular waveguide or coaxial line)
% Modified* Nicolson-Ross-Weir method
% * The transmission medium is air (instead of vacuum) and
%   attenuation due to metal conductivity is taken into account
clear all
close all
clc

% Measured Short and Load (*.s1p)
file_S = '../BM_1.s1p';
file_L =  '../B_1.s1p';
%ref_data = load('../ref_measured.txt');

% MATERIAL SAMPLE PARAMETERS
% Sample cross-section dimensions [m]
SCS = [NaN, 10.03]*1e-3;
% Sample Offset in "Short" and "Load" Cases [m]
L = [0, 0]*1e-3;
% Sample Thickness
Thn = 2.04e-3;

% SAMPLE HOLDER PARAMETERS
% Transmission line type: 'WR' or 'CL'
TL = 'CL';
% Sample holder dimensions [m]
SH = [23.018, 10.03, 6.615]*1e-3;
% Sample holder offset loss at 1 GHz [ohms/s]
OL1GHz = 0e9;

% Band of operation [Hz]
freqs = [3, 20]*1e9;
% Integer number of wavelengths in the sample at the lowest frequency (0, 1, ...)
NRW_root = 0;

% Loss data format: 'tan' for loss tangent
lossf = 'im';
% A degree of a polynomial (eps', eps'', mu', mu'')
PFN = [2, 2, 2, 2];

% Plot properties
FontSize = 18;
LineWidth1 = 1.5;                   % S11, S21, S12, S22
LineWidth2 = 3;                     % fit, air gap correction
MarkerSize = 6;                    % eps, mu
AxisColor = [160, 160, 160]./255;
Color0 = [0, 0, 0];                 % background
Color1 = [.5, .2, .2];              % S11, S21
Color2 = [.2, .2, .5];              % S22, S12
Color3 = [255, 245, 0]./255;        % fit
Color4 = [5, 255, 0]./255;          % air gap correction

[filePath,fileName,fileExt] = fileparts(file_L);
[SS,f] = S_Parameters_Loading(file_S);
SL = S_Parameters_Loading(file_L);
GS = SS(:,1);
GL = SL(:,1);

clear SS SL
I1 = f >= freqs(1);
I2 = f <= freqs(2);
I3 = I1 == I2;
clear I1 I2
f = f(I3);
GS = GS(I3);
GL = GL(I3);
n = length(f);

if 1 == 0
  % TEST PART
  PF1 = polyfit(f, real(GS), 3);
  PF2 = polyfit(f,-imag(GS), 3);
  PF3 = polyfit(f, real(GL), 3);
  PF4 = polyfit(f,-imag(GL), 3);
  GS = polyval(PF1,f)-1i*polyval(PF2,f);
  GL = polyval(PF3,f)-1i*polyval(PF4,f);
  clear PF1, PF2, PF3, PF4;
end


c = 299792458;
mu0 = pi*4e-7;

% The transmission medium (air [2, 3])
epsm = 1.000649*ones(n,1);
mum = ones(n,1);

% Calculation of the attenuation constant based on the offset loss [2]
if strcmp(TL,'WR')
  if SH(1) >= SH(2)
    disp('TE10 mode');
    disp(' ');
    disp('  Sample holder dimensions, mm:');
    disp(['                            width: ', num2str(SH(1)*1000,6)]);
    disp(['                           height: ', num2str(SH(2)*1000,6)]);
    disp(['                           length: ', num2str(SH(3)*1000,6)]);
    ai = 1/(4*SH(1)^2);
    fc = c./(2*SH(1)*sqrt(1.000649));
    alpham = OL1GHz*(sqrt(1.000649)/c^2/mu0).*sqrt(f./fc).*...
    (1+(2*SH(2)/SH(1))*(fc./f).^2)./sqrt(1-(fc./f).^2);
  else
    disp('Error: Sample holder cross-section dimensions are incorrect');
    return
  end
elseif strcmp(TL,'CL')
  if SH(1) > SH(2)
    disp('TEM mode');
    disp(' ');
    disp('  Sample holder dimensions, mm:');
    disp(['                   outer diameter: ', num2str(SH(1)*1000,6)]);
    disp(['                   inner diameter: ', num2str(SH(2)*1000,6)]);
    disp(['                           length: ', num2str(SH(3)*1000,6)]);
    ai = 0;
    alpham = OL1GHz*sqrt(f*1e-9)*pi*1.000649/(c^2*mu0*log(SH(1)/SH(2)));
  else
    disp('Error: Sample holder cross-section dimensions are incorrect');
    return
  end
else
  disp('Error: Transmission line is not defined');
  return
end

l_g = (c./f)./real(sqrt(epsm.*mum-ai*(c./f).^2));

GS = GS.*exp((2i*pi*sqrt(epsm.*mum./((c./f).^2)-ai)+alpham)*2*L(1));
GL = GL.*exp((2i*pi*sqrt(epsm.*mum./((c./f).^2)-ai)+alpham)*2*L(2));

[eps1,mu1] = NRW_RR(f, GS, GL, -ones(n,1), zeros(n,1), ...
ai, Thn, epsm, mum, alpham, NRW_root);

PF1 = polyfit(f, real(eps1),PFN(1));
PF2 = polyfit(f,-imag(eps1),PFN(2));
PF3 = polyfit(f, real(mu1), PFN(3));
PF4 = polyfit(f,-imag(mu1), PFN(4));

eps_fit = polyval(PF1,f)-1i*polyval(PF2,f).*(polyval(PF2,f)>0);
mu_fit = polyval(PF3,f)-1i*polyval(PF4,f).*(polyval(PF4,f)>0);

% Air gap correction
if strcmp(TL,'WR')
  % [4, eqs (2.19), (2.20)]
  dlt = (SH(2)-SCS(2))/SH(2);
  eps_c = eps_fit*(1-dlt)./(1-eps_fit*dlt./epsm);
  mu_c = (mu_fit-dlt*mum)/(1-dlt);
elseif strcmp(TL,'CL')
  % This may be incorrect
  L1 = log(SCS(2)/SH(2))+log(SH(1)/SCS(1));
  L2 = log(SCS(1)/SCS(2));
  L3 = log(SH(1)/SH(2));
  eps_cr = real(eps_fit)*L2./(L3-real(eps_fit)*L1);
  eps_c = eps_cr+1i*(eps_cr.*imag(eps_fit)./real(eps_fit))...
  *L3./(L3-L1*real(eps_fit).*(1+(-imag(eps_fit)./real(eps_fit)).^2));
  mu_c = real(mu_fit)*(L3-L1)/L2+1i*imag(mu_fit)*L3/L2;
  clear eps_cr
else
  return
end

% Scattering parameters
f1 = figure('Units', 'Normalized', 'Position', [.005 .09 .99 .79],...
'Name','S-parameters','Color',Color0);
subplot(2,1,1)
plot(f/1e9,20*log10(abs(GL)),'-','Color',Color2,'LineWidth',LineWidth1);
hold on
plot(f/1e9,20*log10(abs(GS)),'-','Color',Color1,'LineWidth',LineWidth1);
title('\rm\itS\rm_{11}, \itS\rm_{22}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
ylabel('Magnitude, dB','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,1,2)
plot(f/1e9,angle(GL)*180/pi,'-','Color',Color2,'LineWidth',LineWidth1);
hold on
plot(f/1e9,angle(GS)*180/pi,'-','Color',Color1,'LineWidth',LineWidth1);
ylim([-180 180]);
title('\rm\itS\rm_{11}, \itS\rm_{22}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
ylabel('Phase, deg','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'YTick',-180:90:180,'FontSize',FontSize,...
'Color',Color0,'YColor',AxisColor,'XColor',AxisColor);

% Permittivity & permeability
f2 = figure('Units', 'Normalized', 'Position', [.005 .09 .99 .79],...
'Name','Permittivity & permeability','Color',Color0);
subplot(2,2,1)
plot(f/1e9,real(eps1),'.','Color',Color1,'MarkerSize',MarkerSize);
hold on
%plot(ref_data(:,1),ref_data(:,2),'.','Color',Color4,'MarkerSize',MarkerSize);
hold on
%plot(ref_data(:,1),ref_data(:,6),'.','Color',Color4,'MarkerSize',MarkerSize);
hold on
plot(f/1e9,real(eps_fit),'-','Color',Color3,'LineWidth',LineWidth2);
hold on
%plot(f/1e9,real(eps_c),'-','Color',Color4,'LineWidth',LineWidth2);
title('\epsilon_r^\prime','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)
%ylim([8, 18])

subplot(2,2,2)
if isequal(lossf,'tan') == 1
plot(f/1e9,-imag(eps1)./real(eps1),'.','Color',Color1,'MarkerSize',MarkerSize);
hold on
plot(f/1e9,-imag(eps_fit)./real(eps_fit),'-','Color',Color3,'LineWidth',LineWidth2);
hold on
%plot(f/1e9,-imag(eps_c)./real(eps_c),'-','Color',Color4,'LineWidth',LineWidth2);
title('tan\delta_e','FontSize',FontSize,'Color',AxisColor);
else
plot(f/1e9,-imag(eps1),'.','Color',Color1,'MarkerSize',MarkerSize);
hold on
%plot(ref_data(:,1),ref_data(:,3),'.','Color',Color4,'MarkerSize',MarkerSize);
hold on
%plot(ref_data(:,1),ref_data(:,7),'.','Color',Color4,'MarkerSize',MarkerSize);
hold on
plot(f/1e9,-imag(eps_fit),'-','Color',Color3,'LineWidth',LineWidth2);
hold on
%plot(f/1e9,-imag(eps_c),'-','Color',Color4,'LineWidth',LineWidth2);
title('\epsilon_r^{\prime\prime}','FontSize',FontSize,'Color',AxisColor);
end
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)
%ylim([0, 4])

subplot(2,2,3)
plot(f/1e9,real(mu1),'.','Color',Color1,'MarkerSize',MarkerSize);
hold on
%plot(ref_data(:,1),ref_data(:,4),'.','Color',Color4,'MarkerSize',MarkerSize);
hold on
%plot(ref_data(:,1),ref_data(:,8),'.','Color',Color4,'MarkerSize',MarkerSize);
hold on
plot(f/1e9,real(mu_fit),'-','Color',Color3,'LineWidth',LineWidth2);
hold on
%plot(f/1e9,real(mu_c),'-','Color',Color4,'LineWidth',LineWidth2);
title('\mu_r^{\prime}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)
%ylim([.0, 2])

subplot(2,2,4)
if isequal(lossf,'tan') == 1
plot(f/1e9,-imag(mu1)./real(mu1),'.','Color',Color1,'MarkerSize',MarkerSize);
hold on
plot(f/1e9,-imag(mu_fit)./real(mu_fit),'-','Color',Color3,'LineWidth',LineWidth2);
hold on
%plot(f/1e9,-imag(mu_c)./real(mu_c),'-','Color',Color4,'LineWidth',LineWidth2);
title('tan\delta_m','FontSize',FontSize,'Color',AxisColor);
else
plot(f/1e9,-imag(mu1),'.','Color',Color1,'MarkerSize',MarkerSize);
hold on
%plot(ref_data(:,1),ref_data(:,5),'.','Color',Color4,'MarkerSize',MarkerSize);
hold on
%plot(ref_data(:,1),ref_data(:,9),'.','Color',Color4,'MarkerSize',MarkerSize);
hold on
plot(f/1e9,-imag(mu_fit),'-','Color',Color3,'LineWidth',LineWidth2);
hold on
%plot(f/1e9,-imag(mu_c),'-','Color',Color4,'LineWidth',LineWidth2);
title('\mu_r^{\prime\prime}','FontSize',FontSize,'Color',AxisColor);
end
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)
%ylim([-.5, 1])

%return

%dlmwrite('eps_mu_1.txt',...
%[f real(eps1) -imag(eps1) real(mu1) -imag(mu1)],'\t');

dlmwrite(strcat('C:\Documents and Settings\Administrator\Desktop\','S',fileName,'.txt'),...
[f real(eps_fit) -imag(eps_fit) real(mu_fit) -imag(mu_fit)],'\t');

% REFERENCES

% [1]   Measuring the Permittivity and Permeability of Lossy Materials: Solids,
%       Liquids, Metals, Building Materials, and Negative-Index Materials - 
%       NIST Technical Note 1536. J. B. Jarvis, M. D. Janezic, B. F. Riddle,
%       R. T. Johnk, P. Kabos, C. L. Holloway, R. G. Geyer, C. A. Grosvenor.

% [2]   Agilent Technologies. Specifying Calibration Standards and Kits
%       for Agilent Vector Network Analyzers - Application Note 1287-11.

% [3]   Keysight Technologies. Specifying Calibration Standards and Kits
%       for Keysight Vector Network Analyzers - Application Note.

% [4]   Larsson, C., Sjöberg, D., & Elmkvist, L. (2010). Waveguide 
%       measurements of the permittivity and permeability at
%       temperatures up to 1000 C. (Technical Report 
%       LUTEDX/(TEAT-7196)/1-22/(2010); Vol. TEAT-7196).

% Air-filled transmission line (rectangular waveguide or coaxial line)
% Modified* Nicolson-Ross-Weir method
% * The transmission medium is air (instead of vacuum) and
%   attenuation due to metal conductivity is taken into account
clear all
close all
clc

% Measured S-parameters (*.s2p)
file = '../20220330.2.s2p';
% ref_data = load('../eps_mu.txt');

% MATERIAL SAMPLE PARAMETERS
% Sample cross-section dimensions [m]
%SCS = [7.01, 3.015]*1e-3;
SCS = [23.018, 10.02]*1e-3;
% Sample position [m]
L = [0.0, 2.26, NaN]*1e-3;

% SAMPLE HOLDER PARAMETERS
% Transmission line type: 'WR' or 'CL'
TL = 'WR';
% Sample holder dimensions [m]
%SH = [7.02, 3.015, 23.229]*1e-3;
SH = [23.018, 10.03, 6.617]*1e-3;
% Sample holder offset loss at 1 GHz [ohms/s]
%OL1GHz = 1.2e9;
OL1GHz = 0;

% Band of operation [Hz] (used for calculating sample position
% and polynomial curve fitting)
freqs = [8, 12]*1e9;
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
MarkerSize = 9;                    % eps, mu
AxisColor = [160, 160, 160]./255;
Color0 = [0, 0, 0];                 % background
Color1 = [.5, .2, .2];              % S11, S21
Color2 = [.2, .2, .5];              % S22, S12
Color3 = [255, 245, 0]./255;        % fit
Color4 = [5, 255, 0]./255;          % air gap correction

[filePath,fileName,fileExt] = fileparts(file);
[S,f] = S_Parameters_Loading(file);
%S(1601:end,:) = [];
%f(1601:end,:) = [];
%S(1:79,:) = [];
%f(1:79,:) = [];
S11 = S(:,1);
S21 = S(:,2);
S12 = S(:,3);
S22 = S(:,4);
n = length(f);
clear S
I1 = f >= freqs(1);
I2 = f <= freqs(2);
I3 = I1 == I2;
clear I1 I2

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
x3 = median((unwrap(angle(S22(I3)./S11(I3)))).*l_g(I3)./(4*pi));

nanL = isnan(L);
fnanL = find(nanL);
snanL = sum(nanL);

if snanL == 2
  if nanL(2) == 0
    SL2 = ' (user-defined)';
    L(3) = (SH(3)-L(2)-x3)/2;
    SL3 = ' (calculated)';
    L(1) = x3+L(3);
    SL1 = ' (calculated)';
  else
    if nanL(1) == 0
      SL1 = ' (user-defined)';
      L(3) = L(1)-x3;
      SL3 = ' (calculated)';
    else
      SL3 = ' (user-defined)';
      L(1) = L(3)+x3;
      SL1 = ' (calculated)';
    end
    L(2) = SH(3)-L(1)-L(3);
    SL2 = ' (calculated)';
  end
elseif snanL == 1
  if nanL(1) == 1
    SL2 = ' (user-defined)';
    SL3 = ' (user-defined)';
    L(1) = SH(3)-L(3)-L(2);
    SL1 = ' (calculated)';
  elseif nanL(2) == 1
    SL1 = ' (user-defined)';
    SL3 = ' (user-defined)';
    L(2) = SH(3)-L(3)-L(1);
    SL2 = ' (calculated)';
  else
    SL1 = ' (user-defined)';
    SL2 = ' (user-defined)';    
    L(3) = SH(3)-L(2)-L(1);
    SL3 = ' (calculated)';
  end
else
  disp('Error: Sample position is incorrect');
  return
end

S11 = S11.*exp((2i*pi*sqrt(epsm.*mum./((c./f).^2)-ai)+alpham)*2*L(1));
S22 = S22.*exp((2i*pi*sqrt(epsm.*mum./((c./f).^2)-ai)+alpham)*2*L(3));
S21 = S21.*exp((2i*pi*sqrt(epsm.*mum./((c./f).^2)-ai)+alpham)*(L(1)+L(3)));
S12 = S12.*exp((2i*pi*sqrt(epsm.*mum./((c./f).^2)-ai)+alpham)*(L(1)+L(3)));

disp(' ');
disp('  Material sample position, mm:');
disp(['    from 1st ref. plane to sample: ', num2str(L(1)*1000,4), SL1]);
disp(['                 sample thickness: ', num2str(L(2)*1000,4), SL2]);
disp(['    from 2nd ref. plane to sample: ', num2str(L(3)*1000,4), SL3]);
disp(' ');

[eps1,mu1] = NRW(f,S11,S21,ai,L(2),epsm,mum,alpham,NRW_root);
[eps2,mu2] = NRW(f,S22,S12,ai,L(2),epsm,mum,alpham,NRW_root);

PF1 = polyfit(f(I3),real(eps1(I3)+eps2(I3))./2,PFN(1));
PF2 = polyfit(f(I3),-imag(eps1(I3)+eps2(I3))./2,PFN(2));
PF3 = polyfit(f(I3),real(mu1(I3)+mu2(I3))./2,PFN(3));
PF4 = polyfit(f(I3),-imag(mu1(I3)+mu2(I3))./2,PFN(4));

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
hold on
plot(f/1e9,real(eps_c),'-','Color',Color4,'LineWidth',LineWidth2);                      % eps'_c
%hold on
%plot(ref_data(:,1),ref_data(:,2),'.','Color',[0, 1, 1],'MarkerSize',MarkerSize);
title('\epsilon_r^\prime','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)
%ylim([15, 21])

subplot(2,2,2)
if isequal(lossf,'tan') == 1
plot(f/1e9,-imag(eps2)./real(eps2),'.','Color',Color2,'MarkerSize',MarkerSize);         % tan_e_2
hold on
plot(f/1e9,-imag(eps1)./real(eps1),'.','Color',Color1,'MarkerSize',MarkerSize);         % tan_e_1
hold on
plot(f/1e9,-imag(eps_fit)./real(eps_fit),'-','Color',Color3,'LineWidth',LineWidth2);    % tan_e_fit
hold on
plot(f/1e9,-imag(eps_c)./real(eps_c),'-','Color',Color4,'LineWidth',LineWidth2);        % tan_e_c
title('tan\delta_e','FontSize',FontSize,'Color',AxisColor);
else
plot(f/1e9,-imag(eps2),'.','Color',Color2,'MarkerSize',MarkerSize);                     % eps''_2
hold on
plot(f/1e9,-imag(eps1),'.','Color',Color1,'MarkerSize',MarkerSize);                     % eps''_1
hold on
plot(f/1e9,-imag(eps_fit),'-','Color',Color3,'LineWidth',LineWidth2);                   % eps''_fit
hold on
plot(f/1e9,-imag(eps_c),'-','Color',Color4,'LineWidth',LineWidth2);                     % eps''_c
%hold on
%plot(ref_data(:,1),ref_data(:,3),'.','Color',[0, 1, 1],'MarkerSize',MarkerSize);
title('\epsilon_r^{\prime\prime}','FontSize',FontSize,'Color',AxisColor);
end
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)
ylim([-5, 15])

subplot(2,2,3)
plot(f/1e9,real(mu2),'.','Color',Color2,'MarkerSize',MarkerSize);                       % mu'_2
hold on
plot(f/1e9,real(mu1),'.','Color',Color1,'MarkerSize',MarkerSize);                       % mu'_1
hold on
plot(f/1e9,real(mu_fit),'-','Color',Color3,'LineWidth',LineWidth2);                     % mu'_fit
hold on
plot(f/1e9,real(mu_c),'-','Color',Color4,'LineWidth',LineWidth2);                       % mu'_c
%hold on
%plot(ref_data(:,1),ref_data(:,4),'.','Color',[0, 1, 1],'MarkerSize',MarkerSize);
title('\mu_r^{\prime}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)
ylim([-5, 15])

subplot(2,2,4)
if isequal(lossf,'tan') == 1
plot(f/1e9,-imag(mu2)./real(mu2),'.','Color',Color2,'MarkerSize',MarkerSize);           % tan_m_2
hold on
plot(f/1e9,-imag(mu1)./real(mu1),'.','Color',Color1,'MarkerSize',MarkerSize);           % tan_m_1
hold on
plot(f/1e9,-imag(mu_fit)./real(mu_fit),'-','Color',Color3,'LineWidth',LineWidth2);      % tan_m_fit
hold on
plot(f/1e9,-imag(mu_c)./real(mu_c),'-','Color',Color4,'LineWidth',LineWidth2);          % tan_m_c
title('tan\delta_m','FontSize',FontSize,'Color',AxisColor);
else
plot(f/1e9,-imag(mu2),'.','Color',Color2,'MarkerSize',MarkerSize);                      % mu''_2
hold on
plot(f/1e9,-imag(mu1),'.','Color',Color1,'MarkerSize',MarkerSize);                      % mu''_1
hold on
plot(f/1e9,-imag(mu_fit),'-','Color',Color3,'LineWidth',LineWidth2);                    % mu''_fit
hold on
plot(f/1e9,-imag(mu_c),'-','Color',Color4,'LineWidth',LineWidth2);                      % mu''_c
%hold on
%plot(ref_data(:,1),ref_data(:,5),'.','Color',[0, 1, 1],'MarkerSize',MarkerSize);
title('\mu_r^{\prime\prime}','FontSize',FontSize,'Color',AxisColor);
end
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',Color0,...
'YColor',AxisColor,'XColor',AxisColor)
%ylim([0, .8])

%return
MMM = [f*1e-9 real(eps_c) -imag(eps_c) real(mu_c) -imag(mu_c),...
real(eps1) -imag(eps1) real(mu1) -imag(mu1),...
real(eps2) -imag(eps2) real(mu2) -imag(mu2)];
dlmwrite(strcat(filePath,'/',fileName,'_',num2str(L(2)*1000),'_',...
num2str(SH(3)*1000),'_em.txt'),...
MMM,'\t');

%dlmwrite('../ref_measured.txt',...
%[f*1e-9 real(eps1) -imag(eps1) real(mu1) -imag(mu1) ...
%real(eps2) -imag(eps2) real(mu2) -imag(mu2)],'\t');

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

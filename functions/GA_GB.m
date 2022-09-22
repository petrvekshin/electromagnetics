clear all
close all
clc
c = 299792458;
mu0 = pi*4e-7;

f = [1e9:1e7:20e9]';
n = length(f);
%beta = 2i*pi*f*sqrt(1.000649)./c;
ai = 0;
alpha = 0*ones(n,1);
l2 = (c./f).^2;

% Reference media (air)
eps1 = linspace(1.000649,1.000649,n)'-1i*linspace(0,0,n)';
mu1 = linspace(1,1,n)'-1i*linspace(0,0,n)';
Z1 = sqrt(mu1./eps1);

% Material sample
eps2 = linspace(14,17,n)'-1i*linspace(5.5,1.5,n)';
mu2 = linspace(1.8,.6,n)'-1i*linspace(.4,.25,n)';
Z2 = sqrt(mu2./eps2);
d2 = .005;

% Helping slab
eps3 = linspace(3.7,3.6,n)'-1i*linspace(.08,.07,n)';
mu3 = linspace(1,1,n)'-1i*linspace(0,0,n)';
Z3 = sqrt(mu3./eps3);
d3 = .006;

G = (Z2-Z1)./(Z2+Z1);
P = exp(-2i*pi*f.*sqrt(eps2.*mu2)*d2./c);

% Reflection coefficients
GLA = -ones(n,1);                                                       %           Z1 | PEC
GA = Reflection(f, d2, eps2, mu2, Z1, zeros(n,1));                      %      Z1 | Z2 | PEC
GLB = Reflection(f, d3, eps3, mu3, Z1, zeros(n,1));                     %      Z1 | Z3 | PEC
GB = Reflection(f, [d3 d2], [eps3 eps2], [mu3 mu2], Z1, zeros(n,1));    % Z1 | Z2 | Z3 | PEC

%GA2 = (G+((GLA-G)./(1-G.*GLA)).*P.^2)./(1+G.*((GLA-G)./(1-G.*GLA)).*P.^2);
%GB2 = (G+((GLB-G)./(1-G.*GLB)).*P.^2)./(1+G.*((GLB-G)./(1-G.*GLB)).*P.^2);

% Noise
Na = [.0, .0, .0, .0];
Np = [.0, .0, .0, .0];
N1 = (1-Na(1)+2*Na(1)*rand(n,1)).*exp(1i*pi*(-Np(1)+2*Np(1)*rand(n,1))/180);
N2 = (1-Na(2)+2*Na(2)*rand(n,1)).*exp(1i*pi*(-Np(2)+2*Np(2)*rand(n,1))/180);
N3 = (1-Na(3)+2*Na(3)*rand(n,1)).*exp(1i*pi*(-Np(3)+2*Np(3)*rand(n,1))/180);
N4 = (1-Na(4)+2*Na(4)*rand(n,1)).*exp(1i*pi*(-Np(4)+2*Np(4)*rand(n,1))/180);

GLA = GLA.*N1;
GA = GA.*N2;
GLB = GLB.*N3;
GB = GB.*N4;

x1 = (GB - GA - GLB - GA.*GB.*GLB + GLA + GA.*GB.*GLA -...
 GA.*GLB.*GLA + GB.*GLB.*GLA - sqrt(-4*(GB.*GLA-GA.*GLB).^2 +...
 (GA - GB + GLB + GA.*GB.*GLB - GLA - GA.*GB.*GLA + GA.*GLB.*GLA -...
 GB.*GLB.*GLA).^2))./(2*(GB.*GLA-GA.*GLB));
x2 = (GB - GA - GLB - GA.*GB.*GLB + GLA + GA.*GB.*GLA -...
 GA.*GLB.*GLA + GB.*GLB.*GLA + sqrt(-4*(GB.*GLA-GA.*GLB).^2 +...
 (GA - GB + GLB + GA.*GB.*GLB - GLA - GA.*GB.*GLA + GA.*GLB.*GLA -...
 GB.*GLB.*GLA).^2))./(2*(GB.*GLA-GA.*GLB));

Ia = abs(x1)<1;
G1 = x2;
G1(Ia) = x1(Ia);
clear Ia

T2 = (-GA.*G1.*GLA+GA+G1.^2.*GLA-G1)./((GA.*G1-1).*(G1-GLA));

ln1T = log(abs(1./T2))+1i*unwrap(angle(1./T2)+0*pi);
U = -((ln1T./(4*pi*d2)).^2)+ln1T.*alpha/(4*d2*pi^2)-alpha.^2/(2*pi)^2;
mu = mu1.*sqrt(U).*(1+G1)./((1-G1).*sqrt(eps1.*mu1./l2-ai));
eps = l2.*(ai+U)./mu;

disp([max(abs(G-G1)) median(abs(G-G1))])
%disp([f abs(X) abs(x1) abs(x2)])
disp([max(abs(P.^2-T2)) median(abs(P.^2-T2))])
disp([max(abs(eps2-eps)) median(abs(eps2-eps))]);
disp([max(abs(mu2-mu)) median(abs(mu2-mu))]);

%plot(f*1e-9, abs(x1),'r')
%hold on
%plot(f*1e-9, abs(x2),'b')

% Plot properties
FontSize = 20;
AxisColor = [.95 .95 .95];
LineWidth = 1.5;
Line1Color = [1 .2 .2];
Line2Color = [.2 .2 1];
DarkLine = .6; % Color coefficients of dark lines. DarkLineColor = LineColor*DarkLine


f2 = figure('Position',[963 51 948 900],'Name','Permittivity & permeability','Color',[0 0 0]);
subplot(2,2,1)
plot(f/1e9,real(eps),'-','Color',Line2Color,'LineWidth',LineWidth);
hold on
plot(f/1e9,real(eps2),'-','Color',Line1Color,'LineWidth',LineWidth);
title('\epsilon_r^\prime','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,2,2)
plot(f/1e9,-imag(eps),'-','Color',Line2Color,'LineWidth',LineWidth);
hold on
plot(f/1e9,-imag(eps2),'-','Color',Line1Color,'LineWidth',LineWidth);
title('\epsilon_r^{\prime\prime}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,2,3)
plot(f/1e9,real(mu),'-','Color',Line2Color,'LineWidth',LineWidth);
hold on
plot(f/1e9,real(mu2),'-','Color',Line1Color,'LineWidth',LineWidth);
title('\mu_r^{\prime}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)

subplot(2,2,4)
plot(f/1e9,-imag(mu),'-','Color',Line2Color,'LineWidth',LineWidth);
hold on
plot(f/1e9,-imag(mu2),'-','Color',Line1Color,'LineWidth',LineWidth);
title('\mu_r^{\prime\prime}','FontSize',FontSize,'Color',AxisColor);
xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
grid on
set(gca,'FontSize',FontSize,'Color',[0 0 0],...
'YColor',AxisColor,'XColor',AxisColor)
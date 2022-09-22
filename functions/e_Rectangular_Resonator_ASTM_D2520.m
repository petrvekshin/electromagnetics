% Rectangular resonant cavity (described in ASTM Standard D2520)
clear all
close all
clc

% Sample loaded (*.s2p)
SL = '../Sample_285_2.s2p';

% Empty cavity (*.s2p)
EC = '../Empty_285_2.s2p';

% Sample cross-section area [m^2]
SA = .000875^2*pi;

% Cavity dimensions: width, height, length [m]
D = [28.619, 12.6, 90.706]*1e-3;
%D = [71.786, 34, 99.528]*1e-3;

% Plot properties
FontSize = 18;
LineWidth1 = 2.2;                   % empty cavity
LineWidth2 = 2.2;                   % sample loaded

AxisColor = [160, 160, 160]./255;
Color0 = [0, 0, 0];                 % background
Color1 = [.6, .2, .2];              % S21 of empty cavity
Color2 = [.2, .2, .6];              % S12 of empty cavity
Color3 = [.9, .4, .4];              % S21 of sample loaded
Color4 = [.4, .4, .9];              % S12 of sample loaded
MarkerSize = 22;

% Threshold value [dB]
T = 10*log10(.5);

[E,fe] = S_Parameters_Loading(EC);
[L,fl] = S_Parameters_Loading(SL);
E21 = abs(E(:,2));
E12 = abs(E(:,3));
L21 = abs(L(:,2));
L12 = abs(L(:,3));
ne = length(fe);
nl = length(fl);
clear E L

c = 299792458;
epsa = 1.000649;

if ~(D(1)>=D(2))
  return
end

f_max = min([c/D(1)/sqrt(epsa),c/2/D(2)/sqrt(epsa)]);
p_max = floor(2*D(3)*sqrt(epsa*(f_max./c).^2-1./(4*D(1).^2)));
if p_max < 1
  disp('Resonant frequencies are not found');
  return
end

P = [.5:1:p_max+.5]';
P2 = [1:1:p_max]';
N = length(P);
Fr = c*sqrt(((P*.5./D(3)).^2+1/(2*D(1))^2)./epsa);
F2 = c*sqrt(((P2*.5./D(3)).^2+1/(2*D(1))^2)./epsa);

E = sortrows([Fr,ones(N,3);fe,zeros(ne,1),E21,E12],1);
L = sortrows([Fr,ones(N,3);fl,zeros(nl,1),L21,L12],1);

Ie = find(E(:,2));
Il = find(L(:,2));
ME21 = NaN(N-1,4);
ME12 = NaN(N-1,4);
ML21 = NaN(N-1,4);
ML12 = NaN(N-1,4);

%LC = cell(4,1);
%LC{1,1} = 'S21 (empty cavity)';
%LC{2,1} = 'S12 (empty cavity)';
%LC{3,1} = 'S21 (sample loaded)';
%LC{4,1} = 'S12 (sample loaded)';

for p = 1:N-1
  if Ie(p+1)-Ie(p)>3    % at least 3 points
    VE = E(Ie(p)+1:Ie(p+1)-1,[1,3,4]);
    [ME21(p,1),me21] = max(VE(:,2));
    [ME12(p,1),me12] = max(VE(:,3));
    
    ME21(p,2) = VE(me21,1);
    ME12(p,2) = VE(me12,1);
    
    x1 = VE(1,1);
    x2 = VE(end,1);
    figure('Units', 'Normalized', 'Position', [.005 .09 .99 .79],'Name',...
    strcat('TE10',num2str(p,3)),'Color',Color0);
    plot(VE(:,1)*1e-9,20*log10(VE(:,3)),'-','Color',Color2,'LineWidth',LineWidth1);
    hold on
    plot(VE(:,1)*1e-9,20*log10(VE(:,2)),'-','Color',Color1,'LineWidth',LineWidth1);
    %title(strcat('\rmTE_{10',num2str(p,3),'}'),'FontSize',FontSize,'Color',AxisColor);
    xlabel('Frequency, GHz','FontSize',FontSize,'Color',AxisColor);
    ylabel('Magnitude, dB','FontSize',FontSize,'Color',AxisColor);
    grid on
    set(gca,'FontSize',FontSize,'Color',Color0,'YColor',AxisColor,'XColor',AxisColor)

    if (-1)^p<0
%     f2, S21 (empty)      
      for k = me21:length(VE(:,1))
        if VE(k,2)<ME21(p,1)*10^(T/20)
          ME21(p,4) = interp1([VE(k,2), VE(k-1,2)],...
          [VE(k,1), VE(k-1,1)], ME21(p,1)*10^(T/20));
          break
        end
      end
%     f1, S21, (empty)  
      for k1 = 1:me21
        k = me21+1-k1;
        if VE(k,2)<ME21(p,1)*10^(T/20)
          ME21(p,3) = interp1([VE(k,2), VE(k+1,2)],...
          [VE(k,1), VE(k+1,1)], ME21(p,1)*10^(T/20));
          break
        end
      end
%     f2, S12 (empty)      
      for k = me12:length(VE(:,1))
        if VE(k,3)<ME12(p,1)*10^(T/20)
          ME12(p,4) = interp1([VE(k,3), VE(k-1,3)],...
          [VE(k,1), VE(k-1,1)], ME12(p,1)*10^(T/20));
          break
        end
      end      
%     f1, S12, (empty)  
      for k1 = 1:me12
        k = me12+1-k1;
        if VE(k,3)<ME12(p,1)*10^(T/20)
          ME12(p,3) = interp1([VE(k,3), VE(k+1,3)],...
          [VE(k,1), VE(k+1,1)], ME12(p,1)*10^(T/20));
          break
        end
      end
      hold on
      plot([VE(me12,1),ME12(p,3),ME12(p,4)]*1e-9,...
      20*log10([ME12(p,1),ME12(p,1)*10^(T/20),ME12(p,1)*10^(T/20)]),...
      '.','Color',Color2,'MarkerSize',MarkerSize); 
      hold on
      plot([VE(me21,1),ME21(p,3),ME21(p,4)]*1e-9,...
      20*log10([ME21(p,1),ME21(p,1)*10^(T/20),ME21(p,1)*10^(T/20)]),...
      '.','Color',Color1,'MarkerSize',MarkerSize);
      
      if Il(p+1)-Il(p)>3
        VL = L(Il(p)+1:Il(p+1)-1,[1,3,4]);
        [ML21(p,1),ml21] = max(VL(:,2));
        [ML12(p,1),ml12] = max(VL(:,3));
        
        ML21(p,2) = VL(ml21,1);
        ML12(p,2) = VL(ml12,1);
        
        x1 = min([x1, VL(1,1)]);
        x2 = max([x2, VL(end,1)]);
        hold on
        plot(VL(:,1)*1e-9,20*log10(VL(:,3)),'-','Color',Color4,'LineWidth',LineWidth2);
        hold on
        plot(VL(:,1)*1e-9,20*log10(VL(:,2)),'-','Color',Color3,'LineWidth',LineWidth2);
%       f2, S21 (sample)      
        for k = ml21:length(VL(:,1))
          if VL(k,2)<ML21(p,1)*10^(T/20)
            ML21(p,4) = interp1([VL(k,2), VL(k-1,2)],...
            [VL(k,1), VL(k-1,1)], ML21(p,1)*10^(T/20));
            break
          end
        end
%       f1, S21, (sample)  
        for k1 = 1:ml21
          k = ml21+1-k1;
          if VL(k,2)<ML21(p,1)*10^(T/20)
            ML21(p,3) = interp1([VL(k,2), VL(k+1,2)],...
            [VL(k,1), VL(k+1,1)], ML21(p,1)*10^(T/20));
            break
          end
        end
%       f2, S12 (sample)      
        for k = ml12:length(VL(:,1))
          if VL(k,3)<ML12(p,1)*10^(T/20)
            ML12(p,4) = interp1([VL(k,3), VL(k-1,3)],...
            [VL(k,1), VL(k-1,1)], ML12(p,1)*10^(T/20));
            break
          end
        end
%       f1, S12, (sample)  
        for k1 = 1:ml12
          k = ml12+1-k1;
          if VL(k,3)<ML12(p,1)*10^(T/20)
            ML12(p,3) = interp1([VL(k,3), VL(k+1,3)],...
            [VL(k,1), VL(k+1,1)], ML12(p,1)*10^(T/20));
            break
          end
        end        
     

        hold on
        plot([VL(ml12,1),ML12(p,3),ML12(p,4)]*1e-9,...
        20*log10([ML12(p,1),ML12(p,1)*10^(T/20),ML12(p,1)*10^(T/20)]),...
        '.','Color',Color4,'MarkerSize',MarkerSize);
        hold on
        plot([VL(ml21,1),ML21(p,3),ML21(p,4)]*1e-9,...
        20*log10([ML21(p,1),ML21(p,1)*10^(T/20),ML21(p,1)*10^(T/20)]),...
        '.','Color',Color3,'MarkerSize',MarkerSize);
      end 
    end         % if -1^p<0
    xlim([x1,x2]*1e-9)
    ylim manual
    hold on
    plot([F2(p), F2(p)]*1e-9,[-200, 0],'--','Color',AxisColor,'LineWidth',LineWidth1);
  end
end

eps = NaN(N-1,4);
for J = 1:2:N-1
  eps(J,1) = (D(1)*D(3))*(ME21(J,2)-ML21(J,2))/2/SA/ML21(J,2)+1;
  eps(J,2) = (D(1)*D(3))*(ME12(J,2)-ML12(J,2))/2/SA/ML12(J,2)+1;
  eps(J,3) = (D(1)*D(3))/4/SA*(1/(ML21(J,2)/(ML21(J,4)-ML21(J,3)))...
  -1/(ME21(J,2)/(ME21(J,4)-ME21(J,3))));
  eps(J,4) = (D(1)*D(3))/4/SA*(1/(ML12(J,2)/(ML12(J,4)-ML12(J,3)))...
  -1/(ME12(J,2)/(ME12(J,4)-ME12(J,3))));
end
disp(eps)

% REFERENCES

% [1]   Agilent Technologies. Basics of Measuring the Dielectric
%       Properties of Materials - Application Note.

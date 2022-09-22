clear all
close all
clc

EDF = rand(1,1)*exp(2i*pi*rand(1,1));
ERF = rand(1,1)*exp(2i*pi*rand(1,1));
ESF = rand(1,1)*exp(2i*pi*rand(1,1));
ETF = rand(1,1)*exp(2i*pi*rand(1,1));
ELF = rand(1,1)*exp(2i*pi*rand(1,1));
EXF = rand(1,1)*exp(2i*pi*rand(1,1));

EDR = rand(1,1)*exp(2i*pi*rand(1,1));
ERR = rand(1,1)*exp(2i*pi*rand(1,1));
ESR = rand(1,1)*exp(2i*pi*rand(1,1));
ETR = rand(1,1)*exp(2i*pi*rand(1,1));
ELR = rand(1,1)*exp(2i*pi*rand(1,1));
EXR = rand(1,1)*exp(2i*pi*rand(1,1));

S11A = rand(1,1)*exp(2i*pi*rand(1,1));
S21A = rand(1,1)*exp(2i*pi*rand(1,1));
S12A = rand(1,1)*exp(2i*pi*rand(1,1));
S22A = rand(1,1)*exp(2i*pi*rand(1,1));

S11M = EDF+ERF*(S11A+S21A*S12A*ELF/(1-S22A*ELF))/(1-ESF*(S11A+S21A*S12A*ELF/(1-S22A*ELF)));
S22M = EDR+ERR*(S22A+S21A*S12A*ELR/(1-S11A*ELR))/(1-ESR*(S22A+S21A*S12A*ELR/(1-S11A*ELR)));
S21M = EXF+S21A*ETF/((1-S11A*ESF)*(1-S22A*ELF)-ESF*S21A*S12A*ELF);
S12M = EXR+S12A*ETR/((1-S11A*ELR)*(1-S22A*ESR)-ESR*S21A*S12A*ELR);

S11N = (S11M-EDF)/ERF;
S21N = (S21M-EXF)/ETF;
S12N = (S12M-EXR)/ETR;
S22N = (S22M-EDR)/ERR;

S11 = (S11N*(1+S22N*ESR)-ELF*S21N*S12N)/((1+S11N*ESF)*(1+S22N*ESR)-ELF*ELR*S21N*S12N);
S22 = (S22N*(1+S11N*ESF)-ELR*S21N*S12N)/((1+S11N*ESF)*(1+S22N*ESR)-ELF*ELR*S21N*S12N);
S21 = S21N*(1+S22N*(ESR-ELF))/((1+S11N*ESF)*(1+S22N*ESR)-ELF*ELR*S21N*S12N);
S12 = S12N*(1+S11N*(ESF-ELR))/((1+S11N*ESF)*(1+S22N*ESR)-ELF*ELR*S21N*S12N);

disp(abs([S11A-S11]))
disp(abs([S21A-S21]))
disp(abs([S12A-S12]))
disp(abs([S22A-S22]))



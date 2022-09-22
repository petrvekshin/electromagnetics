function [GA] = One_Port_Cal(G1A,G2A,G3A,G1M,G2M,G3M,GM)
% GA is the actual reflection coefficient of the DUT
% G1A, G2A, G3A are the presumed actual reflection coefficients of the three standards
% G1M, G2M, G3M are the measured values of the standards
% GM is the measured value of the DUT

ED = (G1A.*G2A.*G1M.*G3M-G1A.*G3A.*G1M.*G2M-G1A.*G2A.*G2M.*G3M+G2A.*G3A.*G1M.*G2M+G1A.*G3A.*G2M.*G3M-G2A.*G3A.*G1M.*G3M)./...
(G1A.*G2A.*G1M-G1A.*G2A.*G2M-G1A.*G3A.*G1M+G1A.*G3A.*G3M+G2A.*G3A.*G2M-G2A.*G3A.*G3M);

ES = -(G1A.*G2M-G2A.*G1M-G1A.*G3M+G3A.*G1M+G2A.*G3M-G3A.*G2M)./...
(G1A.*G2A.*G1M-G1A.*G2A.*G2M-G1A.*G3A.*G1M+G1A.*G3A.*G3M+G2A.*G3A.*G2M-G2A.*G3A.*G3M);

ER = (G1A.*G1M.*G2M-G1A.*G1M.*G3M-G2A.*G1M.*G2M+G2A.*G2M.*G3M+G3A.*G1M.*G3M-G3A.*G2M.*G3M)./...
(G1A.*G2A.*G1M-G1A.*G2A.*G2M-G1A.*G3A.*G1M+G1A.*G3A.*G3M+G2A.*G3A.*G2M-G2A.*G3A.*G3M)+ED.*ES;

GA = (GM-ED)./(ES.*(GM-ED)+ER);

end
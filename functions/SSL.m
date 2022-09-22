function [S11a,f] = SSL(S,OS,L,DUT,x)
c = 299792458;
 % measured short
[S11Am,f] = S_Parameters_Loading(S);

% measured offset short
S11Bm = S_Parameters_Loading(OS);

% measured load
S11Cm = S_Parameters_Loading(L);

% measured DUT
S11m = S_Parameters_Loading(DUT);

% an actual short
S11Aa = -1;
% an actual offset short
S11Ba = -exp(-4i*pi*f*sqrt(1.000649)*x./c);
ED = S11Cm;

ER = (S11Bm-ED+S11Ba.*S11Bm-S11Ba.*ED)./S11Ba./(1+(S11Bm-ED)./(ED-S11Am));
ES = ER./(ED-S11Am)-1;;

S11a = (S11m-ED)./(ES.*(S11m-ED)+ER);

end
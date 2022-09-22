function [S] = Rpar2Spar(R)
% Converting R-parameters to S-parameters
% R = (R11 R21 R12 R22)
% S = (S11 S21 S12 S22)

S11 = R(:,3)./R(:,4);
S21 = 1./R(:,4);
S12 = (R(:,1).*R(:,4)-R(:,2).*R(:,3))./R(:,4);
S22 = -R(:,2)./R(:,4);

S = [S11 S21 S12 S22];

end
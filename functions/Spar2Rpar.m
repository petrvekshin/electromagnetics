function [R] = Spar2Rpar(S)
% Converting S-parameters to R-parameters
% S = (S11 S21 S12 S22)
% R = (R11 R21 R12 R22)

R11 = (S(:,2).*S(:,3)-S(:,1).*S(:,4))./S(:,2);
R21 = -S(:,4)./S(:,2);
R12 = S(:,1)./S(:,2);
R22 = 1./S(:,2);

R = [R11 R21 R12 R22];

end
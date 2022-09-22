function [RCS] = RCS_Solid_Angle(M,phi_1,phi_2,k)

[n, m] = size(M);
% n is the number of the azimuth points
I = [1:n-1]';

S1 = sind((phi_1*(n-0.5-I)+phi_2*(I-0.5))/(n-1));
clear I
S = abs([S1; sind(phi_2)]-[sind(phi_1); S1]);
clear S1

V1 = reshape(M,n*m,1);
V2 = repmat(S,m,1);
clear S
V2(1:n) = V2(1:n)*k;
V2(n*(m-1)+1:n*m) = V2(n*(m-1)+1:n*m)*k;
V = [V1 V2];
clear V1 V2
V = sortrows(V,1);

S1 = cumsum(V(:,2));
Sv = S1(n*m);
for I = 1:n*m
    if S1(I) >= Sv/2;
        break;
    end
end

RCS = interp1([S1(I-1) S1(I-1)+V(I,2)/2 S1(I)],...
[(V(I-1,1)+V(I,1))/2 V(I,1) (V(I,1)+V(I+1,1))/2],Sv/2);

end
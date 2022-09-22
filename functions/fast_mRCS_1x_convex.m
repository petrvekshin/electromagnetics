function [mRCS] = fast_mRCS_1x_convex(M,k0,lambda)
% Fast mRCS computation of a convex PEC body using PO with numerical integration

n = length(lambda);

% Wavenumber
wn = 2*pi./lambda;

% Finding visible polygons
V1 = M(:,13).*k0(1,1)+M(:,14).*k0(1,2)+M(:,15).*k0(1,3);
V = V1<0;

mRCS = zeros(1,n);
Hf = zeros(1,n);

% H-field computation
for I = 1:n
    Hf(I) = sum(-M(V,16).*(k0(1,1).*M(V,13)...
        +k0(1,2).*M(V,14)+k0(1,3).*M(V,15)).*exp(-2i.*wn(I).*(k0(1,1).*M(V,10)+...
        k0(1,2).*M(V,11)+k0(1,3).*M(V,12))));
end

% mRCS computation
for I = 1:n
    mRCS(I) = (4*pi/lambda(I)^2)*abs(Hf)^2;
end

end
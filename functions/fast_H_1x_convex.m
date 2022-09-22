function [Hf] = fast_H_1x_convex(M,k0,lambda)

n = length(lambda);

wn = 2*pi./lambda;

V1 = M(:,13).*k0(1,1)+M(:,14).*k0(1,2)+M(:,15).*k0(1,3);
V = V1<0;

mRCS = zeros(1,n);
Hf = zeros(1,n);
for I = 1:n
    Hf(I) = sum(-1i*0.5*M(V,16)./pi.*(k0(1,1).*M(V,13)...
        +k0(1,2).*M(V,14)+k0(1,3).*M(V,15)).*exp(-2i.*wn(I).*(k0(1,1).*M(V,10)+...
        k0(1,2).*M(V,11)+k0(1,3).*M(V,12))));
end

for I = 1:n
    mRCS(I) = (4*pi/lambda(I)^2)*abs(Hf)^2;
end

end
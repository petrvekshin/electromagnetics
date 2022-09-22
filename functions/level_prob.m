function [P] = level_prob(XY,x1,x2,probability)

N = length(XY(:,1));

for I1 = 1:N
    if XY(I1,1) >= x1
        break
    end
end
ind1 = I1;

for I2 = 1:N
    if XY(N+1-I2,1) <= x2
        break
    end
end
ind2 = N+1-I2;

if XY(ind1,1) == x1
    if XY(ind2,1) == x2
        M = XY(ind1:ind2,:);
    else
        M = [XY(ind1:ind2,:);[x2,interp1(XY(:,1),...
        XY(:,2),x2)]];
    end
else
    if XY(ind2,1) == x2
        M = [[x1,interp1(XY(:,1),...
        XY(:,2),x1)];XY(ind1:ind2,1:2)];
    else
        M = [[x1,interp1(XY(:,1),...
        XY(:,2),x1)];...
        [x2,interp1(XY(:,1),...
        XY(:,2),x2)]];
    end
end

V = sort(M(:,2));

for I3 = 1:length(V)
    if F_level(M,V(I3)) <= probability
        break;
    end
end
ind_prob = I3;

if probability == 1
    P = V(1,1);
else
    P = interp1([F_level(M,V(ind_prob)) F_level(M,V(ind_prob-1))],...
        [V(ind_prob) V(ind_prob-1)],probability);
end

end


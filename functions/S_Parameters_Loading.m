function [M,f] = S_Parameters_Loading(file)
% M = [S11 S21 S12 S22]

if file(length(file)-1) == '1'
    n = 6;
else
    n = 9;
end
M = importdata(file,' ',n);
S = M.data;
H = M.colheaders;
clear M
if H{1,1} ~= '#'
    S = importdata(file,' ',0);
end    
M = zeros(length(S(:,1)),n-5);
f = S(:,1);

if size(H,2)<4
    H{1,4} = 0;
end
if isequal(H{1,4},'dB') == 1
    M(:,1) = 10.^(S(:,2)/20).*exp(1i*S(:,3)*pi/180);
    if n == 9
        M(:,2) = 10.^(S(:,4)/20).*exp(1i*S(:,5)*pi/180);
        M(:,3) = 10.^(S(:,6)/20).*exp(1i*S(:,7)*pi/180);
        M(:,4) = 10.^(S(:,8)/20).*exp(1i*S(:,9)*pi/180);
    end
elseif isequal(H{1,4},'MA') == 1
    M(:,1) = S(:,2).*exp(1i*S(:,3)*pi/180);
    if n == 9
        M(:,2) = S(:,4).*exp(1i*S(:,5)*pi/180);
        M(:,3) = S(:,6).*exp(1i*S(:,7)*pi/180);
        M(:,4) = S(:,8).*exp(1i*S(:,9)*pi/180);
    end    
else
    M(:,1) = S(:,2)+1i*S(:,3);
    if n == 9
        M(:,2) = S(:,4)+1i*S(:,5);
        M(:,3) = S(:,6)+1i*S(:,7);
        M(:,4) = S(:,8)+1i*S(:,9);
    end
end

end
clear all
close all
clc

fileName = '.../3d_model.mail';

A = importdata(fileName,' ');

N = A(1,2);
P = A(2:A(1,1)+1,:);
T = A(A(1,1)+2:A(1,1)+N+1,:);
clear A

M = zeros(N,16);

M(:,1:3) = P(T(:,1),:)/1000;
M(:,4:6) = P(T(:,2),:)/1000;
M(:,7:9) = P(T(:,3),:)/1000;
M(:,10) = (M(:,1)+M(:,4)+M(:,7))./3;
M(:,11) = (M(:,2)+M(:,5)+M(:,8))./3;
M(:,12) = (M(:,3)+M(:,6)+M(:,9))./3;
A = sqrt((M(:,1).*(M(:,5)-M(:,8))+M(:,4).*(M(:,8)-M(:,2))+M(:,7).*(M(:,2)-M(:,5))).^2+...
    (M(:,2).*(M(:,6)-M(:,9))+M(:,5).*(M(:,9)-M(:,3))+M(:,8).*(M(:,3)-M(:,6))).^2+...
    (M(:,3).*(M(:,4)-M(:,7))+M(:,6).*(M(:,7)-M(:,1))+M(:,9).*(M(:,1)-M(:,4))).^2);
M(:,16) = A./2;
Ux = (M(:,2).*(M(:,6)-M(:,9))+M(:,5).*(M(:,9)-M(:,3))+M(:,8).*(M(:,3)-M(:,6)))./A;
Uy = (M(:,3).*(M(:,4)-M(:,7))+M(:,6).*(M(:,7)-M(:,1))+M(:,9).*(M(:,1)-M(:,4)))./A;
Uz = (M(:,1).*(M(:,5)-M(:,8))+M(:,4).*(M(:,8)-M(:,2))+M(:,7).*(M(:,2)-M(:,5)))./A;
Nx = (M(:,5)-M(:,2)).*(M(:,9)-M(:,6))-(M(:,6)-M(:,3)).*(M(:,8)-M(:,5));
Ny = (M(:,6)-M(:,3)).*(M(:,7)-M(:,4))-(M(:,4)-M(:,1)).*(M(:,9)-M(:,6));
Nz = (M(:,4)-M(:,1)).*(M(:,8)-M(:,5))-(M(:,5)-M(:,2)).*(M(:,7)-M(:,4));
S = sign(Nx.*Ux+Ny.*Uy+Nz.*Uz);
M(:,13) = S.*Ux;
M(:,14) = S.*Uy;
M(:,15) = S.*Uz;

V = M(:,15)>0;
%M(V,13) = -M(V,13);
%M(V,14) = -M(V,14);
%M(V,15) = -M(V,15);

xc = 0.;
yc = 0.;
zc = 0.0;
I = (M(:,10)-xc).*M(:,13)+(M(:,11)-yc).*M(:,14)+(M(:,12)-zc).*M(:,15)<0;
M(I,13) = -M(I,13);
M(I,14) = -M(I,14);
M(I,15) = -M(I,15);

%M(:,13) = -M(:,13);
%M(:,14) = -M(:,14);
%M(:,15) = -M(:,15);

clear A Nx Ny Nz P T Ux Uy Uz N S fileName V I

patch('XData',([M(:,1) M(:,4) M(:,7)])',...
    'YData',([M(:,2) M(:,5) M(:,8)])',...
    'ZData',([M(:,3) M(:,6) M(:,9)])','FaceColor',[0.1 0.92 1]);
axis equal
hold on
quiver3(M(:,10),M(:,11),M(:,12),...
    M(:,13).*M(:,16),M(:,14).*M(:,16),M(:,15).*M(:,16),'r');

xlabel('X')
ylabel('Y')
zlabel('Z')

R1 = sqrt((M(:,1)-M(:,4)).^2+(M(:,2)-M(:,5)).^2+(M(:,3)-M(:,6)).^2);
R2 = sqrt((M(:,1)-M(:,7)).^2+(M(:,2)-M(:,8)).^2+(M(:,3)-M(:,9)).^2);
R3 = sqrt((M(:,7)-M(:,4)).^2+(M(:,8)-M(:,5)).^2+(M(:,9)-M(:,6)).^2);
disp(['The maximum size of a triangle is ' num2str(max([R1; R2; R3]))]);
disp(['The median size of a triangle is ' num2str(median([R1; R2; R3]))]);

disp([min(min([M(:,1) M(:,4) M(:,7)])) min(min([M(:,2) M(:,5) M(:,8)])) ...
min(min([M(:,3) M(:,6) M(:,9)])); max(max([M(:,1) M(:,4) M(:,7)])) ...
max(max([M(:,2) M(:,5) M(:,8)])) max(max([M(:,3) M(:,6) M(:,9)]))]);

disp(['The surface area is ' num2str(sum(M(:,16)))]);

clear R1 R2 R3
%%%

%patch('XData',([M(I,1) M(I,4) M(I,7)])',...
%    'YData',([M(I,2) M(I,5) M(I,8)])',...
%    'ZData',([M(I,3) M(I,6) M(I,9)])','FaceColor',[1 0.4 0.4]);
%axis equal
%hold on
%patch('XData',([M(~I,1) M(~I,4) M(~I,7)])',...
%    'YData',([M(~I,2) M(~I,5) M(~I,8)])',...
%    'ZData',([M(~I,3) M(~I,6) M(~I,9)])','FaceColor',[0.4 0.4 1]);

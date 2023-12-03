clear
close all

%% Generate airfoil
Au = [0.18362 0.33139 0.2805 0.24597];      % upper surface weights
Al = [-0.18362 -0.20378 -0.17535 -0.12035]; % lower surface weights
af = CST_airfoil(Au,Al,51);

%% Wing geometry
Nhalf = 6;
b = 10;
Lambda = 10;
phi = 10;
twist = linspace(0,0,Nhalf+1);
ys = linspace(0,1,Nhalf+1);
chord = 4*sqrt(1-ys.^2);

[vertex,pctrl,cctrl] = geom2gridnew(b,chord,Lambda,phi,twist);

figure
plot3(vertex(:,1),vertex(:,2),vertex(:,3),'o')
hold on
plot3(pctrl(:,1),pctrl(:,2),pctrl(:,3),'x')
daspect([1 1 1])

%% PLLT
N = 2*Nhalf;

alpha = 10;
uinf = repmat([cosd(alpha) 0 sind(alpha)],N,1);
Vinf = 1;

dl = diff(vertex,1,1);
dA = diff(vertex(:,2)).*cctrl; % panel planform area

vij = zeros(N,N,3);
un = zeros(N,3);
ua = zeros(N,3);
for i = 1:N
    r1 = pctrl - vertex(i,:);
    r2 = pctrl - vertex(i+1,:);
    R1 = vecnorm(r1,2,2);
    R2 = vecnorm(r2,2,2);
    vij(:,i,:) = cross(uinf,r2,2)./(R2.*(R2-dot(uinf,r2,2))) + ...
        (R1+R2).*cross(r2,r2,2)./(R1.*R2.*(R1.*R2+dot(r1,r2,2))) - ...
        cross(uinf,r1,2)./(R1.*(R1-dot(uinf,r1,2)));
    vij(i,i,:) = cross(uinf(1,:),r2(i,:))/(R2(i)*(R2(i)-dot(uinf(1,:),r2(i,:)))) - ...
        cross(uinf(1,:),r1(i,:))/(R1(i)*(R1(i)-dot(uinf(1,:),r1(i,:))));
end

%% AE 451 Aeroelasticity 
% Main Ritz-MLLT solver

clc; clear; close all;

%% Input
% Freestream conditions
rho = 0.002377;
qinf = 10;
Vinf = sqrt(2*qinf/rho);
alpha = 5;

% Structural properties
structprop.EI = 8830/144;
structprop.GJ = 13330/144;
eOc = 0.25;

% Airfoil
% Au = [0.18362 0.33139 0.2805 0.24597];      % upper surface weights
Al = [-0.18362 -0.20378 -0.17535 -0.12035]; % lower surface weights
Au = -Al;
geom.af = CST_airfoil(Au,Al,51);

% Wing geometry
Nhalf = 6;                      % number of horseshoe vortices in a semispan
b = 60/12;                         % wingspan
Lambda = -30;                    % c/4 sweep angle (deg)
phi = 0;                        % dihedral angle (deg)
ys = linspace(0,1,Nhalf+1).';     % non-dimensional semispan coordinate of legs
twist = linspace(0,0,Nhalf+1).';  % twist distribution on ys (deg)
chord = 5/12*linspace(1,1,Nhalf+1).';  % chord distribution on ys

%%
[geom.vertex,geom.pctrl,geom.cctrl] = geom2grid(b,chord,Lambda,phi,twist,ys);

figure
plot3(geom.vertex(:,1),geom.vertex(:,2),geom.vertex(:,3),'o')
hold on
plot3(geom.pctrl(:,1),geom.pctrl(:,2),geom.pctrl(:,3),'x')
daspect([1 1 1])

%% LOOP
figure
axw = gca;
hold on

figure
axt = gca;
hold on

N = 2*Nhalf;
dw = zeros(N+1,1);          % initially zero difference in displacement and twist
dtheta = zeros(Nhalf+1,1);

for i = 1:5
    [geom.vertex,geom.pctrl,geom.cctrl] = geom2grid(b,chord,Lambda,phi,twist+dtheta,ys);
    geom.vertex(:,3) = geom.vertex(:,3) + dw;                    % vertically displace the vertices
    geom.pctrl(:,3) = geom.pctrl(:,3) + 0.5*(dw(1:N)+dw(2:N+1)); % displace control points by interpolated value

    [z,t] = MLLT(geom,alpha,Vinf,rho,eOc);

    % Deformations from Ritz
    [Pben,Ptor] = fastritz(2,geom.vertex(N/2+1:N+1,2)/cosd(Lambda),z(N/2+1:N),t(N/2+1:N),Lambda,structprop);

    % Calculate changes to w and theta at vertex nodes
    dw(N/2+1:N+1) = polyval(Pben,geom.vertex(N/2+1:N+1,2)/cosd(Lambda));
    dw(1:N/2) = flipud(dw(N/2+2:N+1));

    dtheta = polyval(Ptor,geom.vertex(N/2+1:N+1,2)/cosd(Lambda));

    plot(axw,ys,dw(N/2+1:N+1))
    plot(axt,ys,dtheta)
end

% Loop Criteria



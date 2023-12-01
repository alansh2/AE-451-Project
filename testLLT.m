clear
close all

%% Generate airfoil
Au = [0.18362 0.33139 0.2805 0.24597];      % upper surface weights
Al = [-0.18362 -0.20378 -0.17535 -0.12035]; % lower surface weights
af = CST_airfoil(Au,Al,51);

%% Wing geometry
b = 10;         % total wingspan
AR = 7;         % aspect ratio
taper = 0.6;    % taper ratio (tip/chord)
LEsweep = 20;   % leading-edge sweep angle (deg)
dih = 0;        % dihedral angle (deg)

Nvorhalf = 5; % number of horseshoe vortices in a half span
twist = linspace(5,-5,Nvorhalf+1); % half-span twist distribution (deg)

%% Create LLT input points from the wing geometry
[vertex,pctrl,cctrl] = geom2grid(b,AR,taper,LEsweep,dih,twist);

N = size(pctrl,1);

dl = diff(vertex,1,1);
dA = diff(vertex(:,2)).*cctrl; % panel planform area

alpha = 5;
uinf = repmat([cosd(alpha) 0 sind(alpha)],N,1);
Vinf = 1;
for i = 1:N
    rij = pctrl(i,:) - vertex;
    r = vecnorm(rij,2,2);
    vij(:,(i-1)*3+(1:3)) = cross(uinf,rij(2:N+1,:),2)./(r(2:N+1).*(r(2:N+1)-dot(uinf,rij(2:N+1,:),2)))+...
        (r(1:N)+r(2:N+1)).*cross(rij(1:N,:),rij(2:N+1,:),2)./(r(1:N).*r(2:N+1).*(r(1:N).*r(2:N+1)+dot(rij(1:N,:),rij(2:N+1,:),2)))-...
        cross(uinf,rij(1:N,:),2)./(r(1:N).*(r(1:N)-dot(uinf,rij(1:N,:),2)));
    vij(i,(i-1)*3+(1:3)) = cross(uinf(1,:),rij(i+1,:),2)./(r(i+1).*(r(i+1)-dot(uinf(1,:),rij(i+1,:),2)))-...
        cross(uinf(1,:),rij(i,:),2)./(r(i).*(r(i)-dot(uinf(1,:),rij(i,:),2)));
    un(i,:) = cross(rij(i+1,:),rij(i,:));
    ua(i,:) = rij(i,:) - dl(i,:)/2;
end
un = un./vecnorm(un,2);
ua = ua./vecnorm(ua,2);

% Iterate G
Gamma = sqrt(1-linspace(-1,1,N).^2);
G = reshape(Gamma/Vinf,1,N);
v = uinf + reshape(G*vij/(4*pi),3,N).';
alf = atan2d(dot(v,un,2),dot(v,ua,2));
for i = N:-1:1
    Cl(i) = Panel2D(af,alf(i));
end

ftest = 2*vecnorm(cross(v,dl./dA,2),2,2).*G.';
%fval = mean((2*vecnorm(cross(v,dl./dA,2),2,2).*G.' - Cl).^2)
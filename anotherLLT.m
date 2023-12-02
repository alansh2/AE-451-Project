clear
close all

set(0,'defaultAxesFontName','Times')

%% Generate airfoil
Au = [0.18362 0.33139 0.2805 0.24597];      % upper surface weights
Al = [-0.18362 -0.20378 -0.17535 -0.12035]; % lower surface weights
af = CST_airfoil(Au,Al,51);

%% Wing geometry
b = 10;         % total wingspan
AR = 20/3;         % aspect ratio
taper = 0.5;    % taper ratio (tip/chord)
LEsweep = 10;   % leading-edge sweep angle (deg)
dih = 10;        % dihedral angle (deg)

Nvorhalf = 25; % number of horseshoe vortices in a half span
twist = linspace(2,-2,Nvorhalf+1); % half-span twist distribution (deg)
% twist = sqrt(1-linspace(0,1,Nvorhalf+1));

%% Create LLT input points from the wing geometry
[vertex,pctrl,cctrl] = geom2grid(b,AR,taper,LEsweep,dih,twist);

N = size(pctrl,1);

alpha = 10;
uinf = [cosd(alpha) 0 sind(alpha)];
Vinf = 1;

dl = diff(vertex,1,1);
dA = diff(vertex(:,2)).*cctrl; % panel planform area

vij = zeros(N,N,3);
un = zeros(N,3);
ua = zeros(N,3);
for j = 1:N
    r2 = pctrl(j,:) - vertex(1,:);
    R2 = vecnorm(r2);
    for i = 1:N
        r1 = r2;
        R1 = R2;
        r2 = pctrl(j,:) - vertex(i+1,:);
        R2 = vecnorm(r2);
        if i == j
            vij(i,j,:) = cross(uinf,r2)/(R2*(R2-dot(uinf,r2)))+...
                (R1+R2)*cross(r1,r2)/(R1*R2*(R1*R2+dot(r1,r2)))-...
                cross(uinf,r1)/(R1*(R1-dot(uinf,r1)));
            un(j,:) = cross(r2,r1);
            ua(j,:) = r1 - dl(j,:)/2;
            % un(j,:) = cross(ua(j,:),dl(j,:));
        else
            vij(i,j,:) = cross(uinf,r2)/(R2*(R2-dot(uinf,r2)))-...
                cross(uinf,r1)/(R1*(R1-dot(uinf,r1)));
        end
    end
end
un = un./vecnorm(un,2); % make direction vectors unit
ua = ua./vecnorm(ua,2);

%% Iterate G
Cl1 = Panel2D(af,0);
Cl2 = Panel2D(af,1);
a0 = (Cl2 - Cl1)*180/pi; % 2D lift curve slope
alfZL = -Cl1/a0;
A = diag(2*vecnorm(cross(repmat(uinf,N,1),dl./dA,2),2,2),0) - a0/(4*pi)*sum(vij.*reshape(un,1,N,3),3).';
RHS = a0*(un*uinf.' - alfZL);
G = A \ RHS;

zeta = dl./dA;
[G,E] = solvegamma(G,alpha,vij,un,ua,zeta,af);

% Check answer for Gamma/Vinf
v = zeros(N,3);
for i = 1:3
    v(:,i) = uinf(i) + vij(:,:,i)*G/(4*pi);
end
alf = atan2d(dot(v,un,2),dot(v,ua,2));
Cl = zeros(N,1);
for i = 1:N
    Cl(i) = Panel2D(af,alf(i));
end
ftest = 2*vecnorm(cross(v,zeta,2),2,2).*G;

S = b^2/AR;
CL = sum(Cl.*cctrl.*dl(:,2))/S;

%% Plot
figure
plot(pctrl(:,2),G)

%%
wmodep = cubicdisplacement(dl(N/2+1:N),ones(N/2,1));
% wmodep = cubicdisplacement(dl(N/2+1:N),[ones(floor(N/4),1);-ones(ceil(N/4),1)]);

figure
plot(pctrl(N/2+1:N,2),polyval(wmodep,pctrl(N/2+1:N,2)))
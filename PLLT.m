clear
close all

%% Generate airfoil
Au = [0.18362 0.33139 0.2805 0.24597];      % upper surface weights
Al = [-0.18362 -0.20378 -0.17535 -0.12035]; % lower surface weights
% Au = -Al;
af = CST_airfoil(Au,Al,51);

%% Wing geometry
Nhalf = 10;                     % number of horseshoe vortices in a semispan
b = 10;                         % wingspan
Lambda = 20;                    % c/4 sweep angle (deg)
phi = 0;                       % dihedral angle (deg)
ys = linspace(0,1,Nhalf+1);     % non-dimensional semispan coordinate of legs
% ys = 0.5*cos(linspace(pi,0,Nhalf+1))+0.5;
twist = linspace(0,0,Nhalf+1);  % twist distribution on ys (deg)
% twist = 0.6*(2/pi^2*sqrt(1-ys.^2)+1/pi/8)*180/pi;
% chord = 2*sqrt(1-ys.^2);         % chord distribution on ys
chord = linspace(2,2,Nhalf+1);

[vertex,pctrl,cctrl] = geom2grid(b,chord,Lambda,phi,twist,ys);

figure
plot3(vertex(:,1),vertex(:,2),vertex(:,3),'o')
hold on
plot3(pctrl(:,1),pctrl(:,2),pctrl(:,3),'x')
daspect([1 1 1])

%% PLLT
N = 2*Nhalf;

alpha = 10;
uinf = [cosd(alpha) 0 sind(alpha)];
Vinf = 1;

dl = diff(vertex,1,1);
dA = diff(vertex(:,2)).*cctrl; % panel planform area
zeta = dl./dA;

vij = zeros(N,N,3);
un = zeros(N,3);
ua = zeros(N,3);
for j = 1:N
    r2 = pctrl(j,:) - vertex(1,:);
    R2 = vecnorm(r2);
    for i = 1:N
        % r1 = pctrl(j,:) - vertex(i,:);
        r1 = r2;
        r2 = pctrl(j,:) - vertex(i+1,:);
        % R1 = vecnorm(r1);
        R1 = R2;
        R2 = vecnorm(r2);
        if i == j
            vij(i,j,:) = cross(uinf,r2)/(R2*(R2-dot(uinf,r2))) +...
                (R1+R2)*cross(r1,r2)/(R1*R2*(R1*R2+dot(r1,r2))) -...
                cross(uinf,r1)/(R1*(R1-dot(uinf,r1)));
            un(i,:) = cross(r2,r1);
            ua(i,:) = r1 - 0.5*dl(i,:);
        else
            vij(i,j,:) = cross(uinf,r2)/(R2*(R2-dot(uinf,r2))) -...
                cross(uinf,r1)/(R1*(R1-dot(uinf,r1)));
        end
    end
end
un = un./vecnorm(un,2,2); % make direction vectors unit
ua = ua./vecnorm(ua,2,2);

%% Iterate G
Cl1 = Panel2D(af,0);
Cl2 = Panel2D(af,1);
a0 = (Cl2 - Cl1)*180/pi; % 2D lift curve slope
alfZL = -Cl1/a0;
A = diag(2*vecnorm(cross(repmat(uinf,N,1),zeta,2),2,2),0) - a0/(4*pi)*sum(vij.*reshape(un,1,N,3),3).';
RHS = a0*(un*uinf.' - alfZL);
G0 = A \ RHS;

[G,E] = solvegamma(G0,alpha,vij,un,ua,zeta,af);

% Check answer for Gamma/Vinf
v = zeros(N,3);
for i = 1:3
    v(:,i) = uinf(i) + vij(:,:,i)*G/(4*pi);
end
alf = atan2d(dot(v,un,2),dot(v,ua,2));
Cl = zeros(N,1);
Cm = zeros(N,1);
for i = 1:N
    [Cl(i),~,Cm(i),~,~] = Panel2D(af,alf(i));
end
ftest = 2*vecnorm(cross(v,zeta,2),2,2).*G;

% S = b^2/AR;
% CL = sum(Cl.*cctrl.*dl(:,2))/S;

%% Plot
figure
plot(pctrl(:,2),G)
hold on
plot(pctrl(:,2),G0,'--')
plot(ys*b/2,G0(Nhalf)*sqrt(1-ys.^2))

figure
plot(pctrl(:,2),Cl.*cctrl/mean(cctrl))
hold on
plot(pctrl(:,2),Cl,'--')

%% Ritz method
figure
hold on

ys = linspace(0,1,51);
for i = [1 2 5 10]
    % Use Cl for now
    P = rrcfpoly(i,vertex(Nhalf+1:N+1,2)/cosd(Lambda),Cl(Nhalf+1:N),Cm(Nhalf+1:N));
    w = zeros(1,51);
    for j = 1:i
        w = w + polyval(flipud(P(1:j+4,j)),ys);
    end
    plot(ys,w,'DisplayName',sprintf('%d modes',i))
end
legend('Location','northwest')
xlabel('y/b')
ylabel('w')
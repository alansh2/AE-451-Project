clear
close all

%% Generate airfoil
Au = [0.18362 0.33139 0.2805 0.24597];      % upper surface weights
Al = [-0.18362 -0.20378 -0.17535 -0.12035]; % lower surface weights
% Au = -Al;
af = CST_airfoil(Au,Al,51);

structprop.EI = 8830/144;
structprop.GJ = 13330/144;
e = 0.25;

qinf = 20;
alpha = 5;

%% Wing geometry
Nhalf = 8;                     % number of horseshoe vortices in a semispan
b = 5;                         % wingspan
Lambda = 30;                    % c/4 sweep angle (deg)
phi = 0;                       % dihedral angle (deg)
ys = linspace(0,1,Nhalf+1).';     % non-dimensional semispan coordinate of legs
% ys = 0.5*cos(linspace(pi,0,Nhalf+1).')+0.5;
twist = linspace(0,0,Nhalf+1).';  % twist distribution on ys (deg)
% twist = 0.6*(2/pi^2*sqrt(1-ys.^2)+1/pi/8)*180/pi;
% chord = 2*sqrt(1-ys.^2);         % chord distribution on ys
chord = 5/12*linspace(1,1,Nhalf+1).';

[vertex,pctrl,cctrl] = geom2grid(b,chord,Lambda,phi,twist,ys);

figure
plot3(vertex(:,1),vertex(:,2),vertex(:,3),'o')
hold on
plot3(pctrl(:,1),pctrl(:,2),pctrl(:,3),'x')
daspect([1 1 1])

%% PLLT
N = 2*Nhalf;

uinf = [cosd(alpha) 0 sind(alpha)];

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
Cl = zeros(Nhalf,1);
Cm = zeros(Nhalf,1);
for i = 0:Nhalf-1
    [Cl(i+1),~,Cm(i+1),~,~] = Panel2D(af,alf(Nhalf-i));
end
ftest = 2*vecnorm(cross(v,zeta,2),2,2).*G;

z = qinf*Cl.*cctrl(Nhalf+1:N);
t = qinf*(Cm + Cl*e).*cctrl(Nhalf+1:N).^2;

%% Ritz method
figure
ax1 = gca;
hold on
xlabel('y/b')
ylabel('w')

figure
ax2 = gca;
hold on
xlabel('y/b')
ylabel('twist')

for i = [1 2 5]
    [Pb,Pt] = fastritz(i,vertex(Nhalf+1:N+1,2)/cosd(Lambda),z,t,Lambda,structprop,'lstsq');
    dw = polyval(Pb,vertex(Nhalf+1:N+1,2)/cosd(Lambda));
    dtheta = polyval(Pt,vertex(Nhalf+1:N+1,2)/cosd(Lambda));
    plot(ax1,ys,dw,'DisplayName',sprintf('%d modes',i))
    plot(ax2,ys,dtheta,'DisplayName',sprintf('%d modes',i))
end
legend(ax1,'Location','northwest')
legend(ax2,'Location','northwest')

%% Update geometry to deformed state
[vertex,pctrl,cctrl] = geom2grid(b,chord,Lambda,phi,twist+dtheta,ys);
vertex(Nhalf+1:N+1,3) = vertex(Nhalf+1:N+1,3) + polyval(Pb,max(0,e*chord*sind(Lambda)+vertex(Nhalf+1:N+1,2)/cosd(Lambda)));
vertex(1:Nhalf,3) = flipud(vertex(Nhalf+2:N+1,3));
pctrl(:,3) = pctrl(:,3) + polyval(Pb,max(0,abs(pctrl(:,2))/cosd(Lambda)-(e*cctrl-(pctrl(:,1)-0.5*(vertex(1:N,1)+vertex(2:N+1,1))))*sind(Lambda)));

figure
plot3(vertex(:,1),vertex(:,2),vertex(:,3),'o')
hold on
plot3(pctrl(:,1),pctrl(:,2),pctrl(:,3),'x')
daspect([1 1 1])
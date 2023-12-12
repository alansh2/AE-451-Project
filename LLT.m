function [z,t] = LLT(geom,qinf,alpha,e)
N = size(geom.pctrl,1);

uinf = [cosd(alpha) 0 sind(alpha)];

dl = diff(geom.vertex,1,1);
dA = diff(geom.vertex(:,2)).*geom.cctrl; % panel planform area
zeta = dl./dA;

vij = zeros(N,N,3);
un = zeros(N,3);
ua = zeros(N,3);
for j = 1:N
    r2 = geom.pctrl(j,:) - geom.vertex(1,:);
    R2 = vecnorm(r2);
    for i = 1:N
        % r1 = geom.pctrl(j,:) - geom.vertex(i,:);
        r1 = r2;
        r2 = geom.pctrl(j,:) - geom.vertex(i+1,:);
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

% Iterate G
Cl1 = Panel2D(geom.af,0);
Cl2 = Panel2D(geom.af,1);
a0 = (Cl2 - Cl1)*180/pi; % 2D lift curve slope
alfZL = -Cl1/a0;
A = diag(2*vecnorm(cross(repmat(uinf,N,1),zeta,2),2,2),0) - a0/(4*pi)*sum(vij.*reshape(un,1,N,3),3).';
RHS = a0*(un*uinf.' - alfZL);
G0 = A \ RHS;

[G,E] = solvegamma(G0,alpha,vij,un,ua,zeta,geom.af);

% Calculate forces and moments
v = zeros(N,3);
for i = 1:3
    v(:,i) = uinf(i) + vij(:,:,i)*G/(4*pi);
end
alf = atan2d(dot(v,un,2),dot(v,ua,2));
Cl = zeros(N/2,1);
Cm = zeros(N/2,1);
for i = 0:N/2-1
    [Cl(i+1),~,Cm(i+1),~,~] = Panel2D(geom.af,alf(N/2-i));
end
% ftest = 2*vecnorm(cross(v,zeta,2),2,2).*G; % Check answer for Gamma/Vinf

z = qinf*Cl.*geom.cctrl(N/2+1:N);
t = qinf*(Cm + Cl*e).*geom.cctrl(N/2+1:N).^2;
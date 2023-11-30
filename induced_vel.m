function fval = induced_vel(alpha,vertex,ctrlpt,c,Gam,af)
c = reshape(c,N,1);
ti = c(2:N)./c(1:N-1);
MACi = 2/3*c(1:N-1)*(1+ti+ti^2)/(1+ti);
uInf = [cosd(alpha) 0 sind(alpha)];
V = 1;
G = Gam./MACi/V;
v = zeros(3*N,N);
for i = 1:N
    for j = 1:N
        rij = ctrlpt(j,:) - vertex(i:i+1,:);
        r = vecnorm(rij,2,2);
        vij = cross(uInf,rij(2,:))/(r(2)*(r(2) - dot(uInf,rij(2,:))))-...
            cross(uInf,rij(1,:))/(r(1)*(r(1) - dot(uInf,rij(1,:))));
        if i ~= j
            vij = vij + (r(1) + r(2))*(cross(rij(1,:),rij(2,:)))/(r(1)*r(2)*(r(1)*r(2)+dot(rij(1,:),rij(2,:))));
        end
        v((i-1)*3+(1:3),j) = vij * MACi(i)/(4*pi);
    end
end
dl = diff(vertex,1,1);
dA = diff(vertex(:,2)) .* (c(1:N-1)+diff(c)/2);
zeta = MACi.*dl./dA;
vi = reshape(v*G,3,N).'+uInf;
alpha = zeros(N,1);
Cl = zeros(N,1);
for i = 1:N
    rij = ctrlpt(i,:) - vertex(i:i+1,:);
    un = cross(rij(2),rij(1));
    un = un/norm(un);
    ua = ctrlpt(i,:) - (vertex(i,:) + (vertex(i+1,:)-vertex(i,:))/2);
    ua = ua/norm(ua);
    alpha(i) = atand(dot(vi,un)/dot(vi,ua));
    Cl(i) = Panel2D(af,alpha(i));
end
fval = mean((2*vecnorm(cross(vi,zeta,2),2,2).*G - Cl).^2);

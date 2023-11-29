function [Cl,Cd,Cp,xc] = Panel2D(af,alpha)
% af: N-by-2 airfoil coords
% alpha: AoA in degrees
alf = alpha*pi/180;

m = size(af,1) - 1; % number of panels

x = af(:,1);
y = af(:,2);

% find panel angles
dx = x(2:m+1) - x(1:m);
dy = y(2:m+1) - y(1:m);
theta = atan2(dy,dx);

% establish control points
xc = 0.5*(x(1:m) + x(2:m+1));
yc = 0.5*(y(1:m) + y(2:m+1));

% convert control point to local panel coords
xt = xc - x(1:m).';
yt = yc - y(1:m).';

% precompute trig expressions as 1-by-m arrays
costh = cos(theta).';
sinth = sin(theta).';

% find control point coords in CS of all panels
xp = xt.*costh + yt.*sinth;
yp = -xt.*sinth + yt.*costh;
x2 = dx.'.*costh + dy.'.*sinth; % 1-by-m
y2 = 0;

% find theta1, theta2, r1, r2
theta1 = atan2(yp,xp);
theta2 = atan2(yp,xp-x2);
dtheta = theta2 - theta1; % precompute
dtheta(logical(eye(m))) = pi;

r1 = sqrt(xp.^2 + yp.^2);
r2 = sqrt((xp-x2).^2 + yp.^2);
ln = log(r2./r1); % precompute

% compute influence coefficients
ap = yp.*ln + xp.*dtheta;
bp = xp.*ln + x2 - yp.*dtheta;
am = -yp.*ln + (x2 - xp).*dtheta;
bm = (x2 - xp).*ln - x2 + yp.*dtheta;
c = 1./(2*pi*x2);

ua = c.*(am.*costh - bm.*sinth);
ub = c.*(ap.*costh - bp.*sinth);
va = c.*(am.*sinth + bm.*costh);
vb = c.*(ap.*sinth + bp.*costh);

u = [ua, zeros(m,1)] + [zeros(m,1), ub];
v = [va, zeros(m,1)] + [zeros(m,1), vb];

A = zeros(m+1);
A(1:m,:) = -u.*sinth.' + v.*costh.';
A(m+1,1) = 1; % Kutta condition
A(m+1,m+1) = 1;

% compute RHS
b = [cos(alf)*sinth - sin(alf)*costh, 0].';

gamma = A \ b;

B = u.*costh.' + v.*sinth.';

Qt = B*gamma + cos(theta - alf);
Cp = 1 - Qt.^2;
Cn = -Cp.'*dx;
Ca = Cp.'*dy;
Cl = Cn*cos(alf) - Ca*sin(alf);
Cd = Cn*sin(alf) + Ca*cos(alf);

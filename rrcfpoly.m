function P = rrcfpoly(N,ynode,z,t)
% Inputs:
%       N: number of modes (first mode is quadratic)
%   ynode: distribution of nodes that discretizes the beam span
%       z: dimensional lift force distribution (defined at control pts)
%       t: dimensional applied torsion distribution (defined at control pts)
% Polynomials are functions of y/b
% Shape functions are invariant of beam length
% The minimum polynomial order is the order of the GDE
n = length(z);
z = reshape(z,1,n);
t = reshape(t,1,n);
dy = reshape(diff(ynode),n,1);
y = reshape(ynode(1:n),1,n) + 0.5*dy.'; % trapezoidal evaluation points

E = 1;
I = 1;

pd = N + 3; % maximum polynomial degree

A = zeros(pd+1,N);
for i = 1:N
    R(1,:) = (i+1:i+3).*(i:i+2);
    R(2,:) = R(1,:).*(i-1:i+1);
    R = rref(R);
    A(i+2:i+4,i) = [-R(:,3);1];
end

B = A(3:pd+1,:).*(2:pd).'.*(1:pd-1).'; % second derivative

ypow = (y/ynode(n+1)).^((2:pd).'); % precompute necessary powers of y/b

C = zeros(N); % influence coefficients
RHS = zeros(N,1);
for i = 1:N
    Q = spdiags(B(i+2:-1:1,i).'.*ones(pd-1,1),-1-i:0,i+pd,pd-1);
    for j = i:N
        % pmult = Q(1:j+i+3,1:j+2)*B(1:j+2,j); % general polynomial multiplication
        % pmult = [0;pmult./(1:j+i+3).'] % integrate
        pmult = [0;Q(1:j+i+3,1:j+2)*B(1:j+2,j)./(1:j+i+3).']; % one step
        C(i,j) = polyval(flipud(pmult),1); % evaluate at b
    end
    % RHS = -(z.*polyval(flipud(A(1:i+4,i)),y/b)) * dy;
    RHS(i) = -(z.*(A(i+2:i+4,i).'*ypow(i:i+2,:))) * dy; % reduced number of computations
end
C = C + triu(C,1).'; % fill skipped duplicate computations

Ac = (E*I*C) \ RHS;
P = Ac.'.*A;

% Use governing equations to scale coefficients
D = E*I*sparse(P(5:pd+1,:).*((4:pd).*(3:pd-1).*(2:pd-2).*(1:pd-3)).'); % fourth derivative
F = [ones(n,1) (y/ynode(n+1)).'.^(1:pd-4)]*D;
K = (F'*F) \ (F'*(z.')); % solve minimization problem
P = K.'.*P;
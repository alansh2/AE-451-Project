function [Pben,Ptor] = fastritz(N,ynode,z,t,Lam,structprop,method)
% Inputs:
%          N: number of modes (first mode is quadratic)
%      ynode: distribution of nodes that discretizes the beam span
%          z: dimensional lift force distribution (defined at control pts)
%          t: dimensional applied torsion distribution (defined at control pts)
%        Lam: sweep angle of the elastic axis
% structprop: struct with fields {E,I,G,J} for the material and C/S properties
%     method: (optional)
%      - ritz (default): finds coefficients that minimize the functional
%      - lstsq         : finds the least squares solution to the governing equations
% Polynomials are functions of y/b
% Shape functions are invariant of beam length
% The minimum polynomial order is the order of the GDE
% Solutions are stable up to 10 modes
if nargin == 6
    method = 'ritz';
end
if ~ismember(method,{'ritz','lstsq'})
    error('Requested method is not implemented.')
end

n = length(z);
z = reshape(z,n,1);
t = reshape(t,n,1);
dy = reshape(diff(ynode),n,1);
y = reshape(ynode(1:n),n,1) + 0.5*dy; % trapezoidal evaluation points

% Finite difference matrix for first derivative
V1 = ones(n-2,1)./(y(3:n)-y(1:n-2));
D = spdiags([-[V1;1;1] zeros(n,1) [1;1;V1]],-1:1,n,n);
D(1,1:2) = [-1 1]/(y(2)-y(1));
D(n,n-1:n) = [-1 1]/(y(n)-y(n-1));

% Work on bending %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pd = N + 3;

ypow = (y/ynode(n+1)).^(0:pd); % precompute necessary powers of y/b

% The second derivative coefficients follow a nice pattern
B(3,:) = cumsum(2*(3:pd-1)) + 6;
B(2,:) = -2*B(3,:);
B(1,:) = B(3,:);

% Integrate second derivatives to get original polynomials rather than use rref
Pben = [B(1,:)./(1:pd-3);B(2,:)./(3:pd-1);ones(1,N)];
Pben(1:2,:) = Pben(1:2,:)./(2:pd-2);

% Solve for polynomial coefficients
if strcmp(method,'ritz')
    nd = max([1 floor(N/2)]);
    for i = N:-1:1
        Q = spdiags(repmat(B(:,i),1,3).',-2:0,5,3); % multiplication target matrix
        C(i,i:N) = sum(Q*B(:,i:N)./(repmat(2*i-1:N+i-1,5,1)+(0:4).'),1); % evaluating 1 is a simple sum
        C(i+1:N,i) = C(i,i+1:N); % matrix is symmetric
    end
    RHS = (ypow(:,3:pd+1)*spdiags(Pben.',0:-1:-2,pd-1,N).*dy).'*z;
    Ac = (structprop.E*structprop.I*C) \ RHS;
    Pben = flipud(spdiags(Pben.',[-2 -3 -4],nd+4,nd)*Ac(1:nd));
else
    % Use governing equations to find best linear combination of modes
    if N == 1
        M = structprop.E*structprop.I*Pben(3)*24;
    else
        M = structprop.E*structprop.I*spdiags(([B(1,3:N) 0 0;B(2,2:N) 0;B(3,:)].*(2:pd-2).*(1:pd-3)).',-2:0,N,N).';
    end
    F = ypow(:,1:N)*M; % evaluate all polynomials
    K = (F'*F) \ (F'*(z*cosd(Lam)+D*t*sind(Lam)*cosd(Lam))); % solve minimization problem
    Pben = flipud(spdiags(Pben.',[-2 -3 -4],pd+1,N)*K); % adjust magnitudes
end

% Work on twist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pd = N + 1;

% Both the mode and its first derivative have patterned coefficients
B = [-(2:pd);2:pd];
Ptor = [B(1,:)./(1:pd-1);ones(1,N)];

if strcmp(method,'ritz')
    for i = N:-1:1
        % C(i,i:N) = sum([[B(:,i);0] [0;B(:,i)]]*B(:,i:N)./(repmat(2*i-1:N+i-1,3,1)+(0:2).'),1);
        C(i,i:N) = sum((i+1)*[1;-2;1].*(i+1:N+1)./(repmat(2*i-1:N+i-1,3,1)+(0:2).'),1);
        C(i+1:N,i) = C(i,i+1:N);
    end
    RHS = (ypow(:,2:pd+1)*spdiags(Ptor.',[0 -1],pd,N).*dy).'*t;
    Ac = (structprop.G*structprop.J*C) \ RHS;
    Ptor = flipud(spdiags(Ptor.',[-1 -2],nd+2,nd)*Ac(1:nd));
else
    M = structprop.G*structprop.J*spdiags(([B(1,2:N) 0;B(2,:)].*(1:pd-1)).',-1:0,N,N).';
    F = ypow(:,1:N)*M;
    K = (F'*F) \ (F'*(-t*cosd(Lam)^2));
    Ptor = flipud(spdiags(Ptor.',[-1 -2],pd+1,N)*K);
end
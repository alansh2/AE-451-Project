function [coords] = CST_airfoil(Au,Al,Nx,N1,N2)
% Input:
%     Au: upper surface weights
%     Al: lower surface weights
%     Nx: number of points to discretize the chord (output contains 2*Nx-1 points)
%     N1: (optional) first shape parameter
%     N2: (optional) second shape parameter
%
% Output:
% coords: (2*Nx-1)-by-2 array of airfoil coordinates with x in the first column and y in the second column

if nargin == 3
    N1 = 0.5;
    N2 = 1;
end

n = length(Au)-1; % set Bernstein polynomial order to no. weights - 1

% generate x/c grid from TE to LE
x = 0.5*(1 + cos(linspace(0,pi,Nx)))';
% class function
C = x.^N1 .* (1-x).^N2;
% component shape functions
i = 0:n;
K = factorial(n)./factorial(i)./factorial(n-i);
S = K.*x.^i .* (1-x).^(n-i);
% overall shape function
Su = S*Au(:);
Sl = S*Al(:);

yu = C.*Su;
yl = C.*Sl;

coords = [x(1:end-1) yl(1:end-1);flip(x) flip(yu)];

% Kulfan, B. M., "Universal Parametric Geometry Representation Method," Journal of Aircraft, Vol. 45, No. 1, 2008, pp. 142-158.
% https://doi.org/10.2514/1.29958

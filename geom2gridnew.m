function [vertex,pctrl,cctrl] = geom2gridnew(b,c,sweep,dih,twist,varargin)
% Inputs:
% (All spatial distributions are defined along the c/4 line)
%       b: wingspan
%       c: semispan chord distribution
%   sweep: c/4 sweep angle (deg)
%     dih: dihedral angle (deg)
%   twist: semispan twist distribution (deg)
%       y: (optional) semispan distribution of horseshoe legs
% Outputs:
%  vertex: [x y z] coordinates of the horseshoe corners
%   pctrl: [x y z] coordinates of the control points at each panel 3c/4
%   cctrl: chord lengths from sections at control points
N = length(c);
if nargin == 6
    y = reshape(varargin{1},N,1);
    ys = y/y(N);
else
    ys = linspace(0,1,N).';
end
tctrl = reshape(twist(1:N-1) + diff(twist),N-1,1);

vertex = zeros(2*N-1,3);
pctrl = zeros(2*N-2,3);
cctrl = zeros(2*N-2,1);

vertex(N:2*N-1,:) = [tand(sweep)*b/2*ys/cosd(dih) b/2*ys/cosd(dih) zeros(N,1)];

cctrl(N:2*N-2) = c(1:N-1) + 0.5*diff(c);
cctrl(1:N-1) = flipud(cctrl(N:2*N-2));

pctrl(N:2*N-2,:) = [vertex(N:2*N-2,1)+diff(vertex(N:2*N-1,1))+0.5*cctrl(N:2*N-2).*cosd(tctrl) vertex(N:2*N-2,2)+0.5*diff(vertex(N:2*N-1,2)) -0.5*cctrl(N:2*N-2).*sind(tctrl)];

vertex(N:2*N-1,[2 3]) = vertex(N:2*N-1,[2 3])*[cosd(dih) sind(dih);-sind(dih) cosd(dih)];
pctrl(N:2*N-2,[2 3]) = pctrl(N:2*N-2,[2 3])*[cosd(dih) sind(dih);-sind(dih) cosd(dih)];

vertex(1:N-1,:) = flipud(vertex(N+1:2*N-1,:));
vertex(1:N-1,2) = -vertex(1:N-1,2); % y coordinate flips sign
pctrl(1:N-1,:) = flipud(pctrl(N:2*N-2,:));
pctrl(N:2*N-2,2) = -pctrl(N:2*N-2,2);
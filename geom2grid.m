function [vertex,pctrl,cctrl] = geom2grid(b,c,sweep,dih,twist,varargin)
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

vertex = [tand(sweep)*b/2*ys/cosd(dih) b/2*ys/cosd(dih) zeros(N,1)];
cctrl = reshape(c(1:N-1) + 0.5*diff(c),N-1,1);
pctrl = [vertex(1:N-1,1)+0.5*diff(vertex(:,1))+0.5*cctrl.*cosd(tctrl) vertex(1:N-1,2)+0.5*diff(vertex(:,2)) -0.5*cctrl.*sind(tctrl)];

vertex(:,[2 3]) = vertex(:,[2 3])*[cosd(dih) sind(dih);-sind(dih) cosd(dih)];
pctrl(:,[2 3]) = pctrl(:,[2 3])*[cosd(dih) sind(dih);-sind(dih) cosd(dih)];

vertex = [flipud(vertex);vertex(2:N,:)];
vertex(1:N-1,2) = -vertex(1:N-1,2); % y coordinate flips sign
pctrl = [flipud(pctrl);pctrl];
pctrl(1:N-1,2) = -pctrl(1:N-1,2);
cctrl = [flipud(cctrl);cctrl];
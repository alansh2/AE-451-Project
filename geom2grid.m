function [p25,p75,c75out] = geom2grid(b,AR,taper,LEsweep,dih,t25,varargin)
% Inputs:
%       b: span
%      AR: aspect ratio
%   taper: taper ratio (ctip/croot)
% LEsweep: sweep angle of the LE (deg)
%     dih: dihedral angle (deg)
%     t25: twist distribution (deg) sampled at the c/4 span nodes (y25)
%     y25: (optional) non-dimensional span distribution of horseshoe legs
% Outputs:
%     p25: [x y z] coordinates of the horseshoe corners
%     p75: [x y z] coordinates of the control points at each panel 3c/4

Nhalf = length(t25);
t25 = reshape(t25,1,Nhalf);
if nargin == 7
    y25 = reshape(varargin{1},1,Nhalf);
else
    y25 = linspace(0,1,Nhalf); % non-dimensional span coordinates of nodes along c/4 line
end
y75 = y25(1:Nhalf-1) + diff(y25)/2; % for control points along 3c/4 line

croot = b/AR * 2/(1+taper); % trapezoidal planform
ctip = croot*taper;
c75 = croot + y75*(ctip-croot);
c75out = [c75(Nhalf-1:-1:1) c75].';

xLEtip = tand(LEsweep)*b/2;
x25 = croot/4 + y25*(xLEtip+ctip/4-croot/4);
x75 = x25(1:Nhalf-1) + diff(x25)/2 + c75/2;

t75 = interp1(y25,t25,y75,'linear');
z75 = -x75.*sind(t75);
x75 = x75.*cosd(t75);

z25 = linspace(0,b/2*tand(dih),Nhalf);
y75 = y75/cosd(dih); % pre-rotation correction due to dihedral
yz75 = [cosd(dih) -sind(dih);sind(dih) cosd(dih)]*[b/2*y75;z75];

p25 = [x25(Nhalf:-1:2) x25;b/2*[-y25(Nhalf:-1:2) y25];z25(Nhalf:-1:2) z25].';
p75 = [x75(Nhalf-1:-1:1) x75;-yz75(1,Nhalf-1:-1:1) yz75(1,:);yz75(2,Nhalf-1:-1:1) yz75(2,:)].';
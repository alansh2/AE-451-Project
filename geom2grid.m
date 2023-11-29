function [p25,p75] = geom2grid(b,AR,taper,LEsweep,t25,varargin)
% Inputs:
%       b: span
%      AR: aspect ratio
%   taper: taper ratio (ctip/croot)
% LEsweep: sweep angle of the LE (deg)
%     t25: twist distribution (deg) sampled at the c/4 span nodes (y25)
%     y25: (optional) non-dimensional span distribution of horseshoe legs
% Outputs:
%     p25: [x y z] coordinates of the horseshoe corners
%     p75: [x y z] coordinates of the control points at each panel 3c/4

Nhalf = length(t25);
t25 = reshape(t25,1,Nhalf);
if nargin == 6
    y25 = reshape(varargin{1},1,Nhalf);
else
    y25 = linspace(0,1,Nhalf); % non-dimensional span coordinates of nodes along c/4 line
end
y75 = y25(1:Nhalf-1) + diff(y25)/2; % for control points along 3c/4 line

croot = b/AR * 2/(1+taper); % trapezoidal planform
ctip = croot*taper;
c75 = croot + y75*(ctip-croot);

xLEtip = tand(LEsweep)*b/2;
x25 = croot/4 + y25*(xLEtip+ctip/4-croot/4);
x75 = x25(1:Nhalf-1) + diff(x25)/2 + c75/2;

t75 = interp1(y25,t25,y75,'linear');
z75 = -x75.*sind(t75);
x75 = x75.*cosd(t75);

p25 = [[x25(Nhalf:-1:2) x25];b/2*[-y25(Nhalf:-1:2) y25];zeros(1,2*Nhalf-1)].';
p75 = [[x75(Nhalf-1:-1:1) x75];b/2*[-y75(Nhalf-1:-1:1) y75];[z75(Nhalf-1:-1:1) z75]].';
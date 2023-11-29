% AE 451 AEROELASTICITY
% Modern LLT

clc; clear; close all;

%% Wing Geometry

RA = 10;            % Aspect Ratio
c = 1;              % Chord
b = RA*c;           % Wingspan (rectangular)

L = 10;              % Sweep angle [deg]
Phi = 0;            % Dihedral angle [deg]

N = 12;             % Number of horseshoe vortices, control pts (even #)

nodes = [abs(linspace(-b/2,b/2,N+1).*sind(L))+c/4;      % Horseshoe Nodal Points at c/4
    linspace(-b/2,b/2,N+1).*cosd(L)];

P = [abs(linspace(nodes(2,1) + 0.5*diff(nodes(2,1:2)), nodes(2,end) - 0.5*diff(nodes(2,1:2)),N)).*sind(L) + 0.75*c;         % Control points at 3c/4
    linspace(nodes(2,1) + 0.5*diff(nodes(2,1:2)), nodes(2,end) - 0.5*diff(nodes(2,1:2)),N).*cosd(L)];


figure
plot(nodes(2,:), nodes(1,:),'o')
hold on
plot(P(2,:), P(1,:),'o')
hold off
xlabel('y'); ylabel('x');
set(gca, 'YDir','reverse'); set(gca, 'XDir','reverse')
axis equal

%% Induced velocity
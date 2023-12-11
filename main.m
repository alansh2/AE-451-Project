%% AE 451 Aeroelasticity 
% Main Ritz-MLLT solver

clc; clear; close all;

%% Initialize

    % Freestream conditions
    rho = 1;
    Vinf = 1;
    alpha = 5;

    % Structural properties
    E = 1;
    I = 1;
    G = 1;
    J = 1;
    structprop = {E,I,G,J};



%% Generate Initial Geometry

    % Generate airfoil
    Au = [0.18362 0.33139 0.2805 0.24597];      % upper surface weights
    Al = [-0.18362 -0.20378 -0.17535 -0.12035]; % lower surface weights
    % Au = -Al;
    af = CST_airfoil(Au,Al,51);

    Nhalf = 6;                     % number of horseshoe vortices in a semispan
    b = 10;                         % wingspan
    Lambda = 20;                    % c/4 sweep angle (deg)
    phi = 0;                       % dihedral angle (deg)
    ys = linspace(0,1,Nhalf+1);     % non-dimensional semispan coordinate of legs
    % ys = 0.5*cos(linspace(pi,0,Nhalf+1))+0.5;
    twist1 = linspace(0,0,Nhalf+1);  % twist distribution on ys (deg)
    % twist = 0.6*(2/pi^2*sqrt(1-ys.^2)+1/pi/8)*180/pi;
    % chord = 2*sqrt(1-ys.^2);         % chord distribution on ys
    chord = linspace(2,2,Nhalf+1);
    
    [vertex1,pctrl1,cctrl] = geom2grid(b,chord,Lambda,phi,twist1,ys);
    
    % Plot initial geometry
    figure
    plot3(vertex1(:,1),vertex1(:,2),vertex1(:,3),'o')
    hold on
    plot3(pctrl1(:,1),pctrl1(:,2),pctrl1(:,3),'x')
    daspect([1 1 1])

% Initialize w,theta
    % w = zeros(Nhalf,1);
    w = transpose(linspace(0,4,Nhalf+1));
    % theta = zeros(1,Nhalf+1);
    theta = linspace(-20,20,Nhalf+1);
    vertex = vertex1;
    pctrl = pctrl1;

%% LOOP

% Adjust geometry - add onto initial geom!
    % Add bending displacement to vertex nodes
    vertex(Nhalf+2:end,3) = vertex1(Nhalf+2:end,3) + w(2:end);
    vertex(1:Nhalf+1,3) = vertex1(1:Nhalf+1,3) + flip(w);

    % Add twist, define control points
    twist = twist1+ theta;
    tctrl = transpose(twist(1:end-1) +0.5*diff(twist));
    pctrl = [vertex(1:end-1,1)+0.5*diff(vertex(:,1))+0.5*cctrl.*[flip(cosd(tctrl)); cosd(tctrl)] vertex(1:end-1,2)+0.5*diff(vertex(:,2)) vertex(1:end-1,3) + 0.5*diff(vertex(:,3))-0.5*cctrl.*[flip(sind(tctrl)); sind(tctrl)]];


    figure
    plot3(vertex(:,1),vertex(:,2),vertex(:,3),'o')
    hold on
    plot3(pctrl(:,1),pctrl(:,2),pctrl(:,3),'x')
    daspect([1 1 1])



% Aero calculations from MLLT

[z,t] = MLLT(vertex,pctrl,cctrl,alpha,Vinf,rho,af);


% Deflection calculations from Ritz

[Pben,Ptor] = fastritz(N,ynode,z,t,Lambda,structprop,method);
% Calculate w,theta at vertex nodes
w = polyval(Pben,vertex(:,2));
theta = polyval(Ptor,vertex(:,2));


% Loop Criteria



clear
close all

%% Inputs
method = 'lstsq';

structprop.EI = 8830/144;
structprop.GJ = 13330/144;
e = 0.25;

qinf = 20;
alpha = 5;

Nhalf = 8;                              % number of horseshoe vortices in a semispan
b = 5;                                  % wingspan
Lambda = 30;                            % c/4 sweep angle (deg)
phi = 0;                                % dihedral angle (deg)
ys = linspace(0,1,Nhalf+1).';           % non-dimensional semispan coordinate of legs
twist = linspace(0,0,Nhalf+1).';        % twist distribution on ys (deg)
chord = 5/12*linspace(1,1,Nhalf+1).';   % chord distribution on ys

Au = [0.18362 0.33139 0.2805 0.24597];      % upper surface weights
Al = [-0.18362 -0.20378 -0.17535 -0.12035]; % lower surface weights
g0.af = CST_airfoil(Au,Al,51);

%% Process initial geometry
N = 2*Nhalf;
[g0.vertex,g0.pctrl,g0.cctrl] = geom2grid(b,chord,Lambda,phi,twist,ys);

% Align initial geometry with the elastic axis
gb.vertex = (g0.vertex(Nhalf+1:N+1,1:2)-[e*chord(1) 0])*[cosd(Lambda) sind(Lambda);-sind(Lambda) cosd(Lambda)];
gb.pctrl = (g0.pctrl(Nhalf+1:N,1:2)-[e*chord(1) 0])*[cosd(Lambda) sind(Lambda);-sind(Lambda) cosd(Lambda)];
i1 = find(gb.vertex(:,2)>0,1):Nhalf+1; % query positive positions along the span
i2 = find(gb.pctrl(:,2)>0,1):Nhalf;

%% Loop
figure
ax1 = gca;
hold on
xlabel('y/b')
ylabel('w')

figure
ax2 = gca;
hold on
xlabel('y/b')
ylabel('twist')

maxiter = 50;   % maximum allowed iterations before terminating
iter = 0;       % iteration counter
geom = g0;      % set iterable to initial state
xLast = [0;0];  % storage for state of previous iteration 

% Perform iteration 0
[z,t] = LLT(geom,qinf,alpha,e);
[Pben,Ptor] = fastritz(5,geom.vertex(Nhalf+1:N+1,2)/cosd(Lambda),z,t,Lambda,structprop,method);
x = [polyval(Pben,b/2/cosd(Lambda));polyval(Ptor,b/2/cosd(Lambda))];

while (abs(x(1)-xLast(1)) > 1e-3) || (abs(x(2)-xLast(2)) > 1e-3)
    if iter >= maxiter
        warning('Exceeded maximum iterations. Solution did not converge within tolerance.')
        break
    end
    iter = iter + 1;
    xLast = x;

    % Move vertices
    dw = polyval(Pben,gb.vertex(i1,2));
    dtheta = polyval(Ptor,gb.vertex(i1,2));
    geom.vertex(Nhalf+i1,1:2) = [gb.vertex(i1,1).*cosd(dtheta)+g0.vertex(Nhalf+i1,3).*sind(dtheta),...
        gb.vertex(i1,2)]*[cosd(Lambda) -sind(Lambda);sind(Lambda) cosd(Lambda)] + [e*chord(1) 0];
    geom.vertex(Nhalf+i1,3) = -gb.vertex(i1,1).*sind(dtheta) + g0.vertex(Nhalf+i1,3).*cosd(dtheta) + dw;
    geom.vertex(1:Nhalf,3) = flipud(geom.vertex(Nhalf+2:N+1,3));

    plot(ax1,ys,[zeros(i1(1)-1,1);dw])
    plot(ax2,ys,[zeros(i1(1)-1,1);dtheta])

    % Move control points
    dw = polyval(Pben,gb.pctrl(i2,2));
    dtheta = polyval(Ptor,gb.pctrl(i2,2));
    geom.pctrl(Nhalf+i2,1:2) = [gb.pctrl(i2,1).*cosd(dtheta)+g0.pctrl(Nhalf+i2,3).*sind(dtheta),...
        gb.pctrl(i2,2)]*[cosd(Lambda) -sind(Lambda);sind(Lambda) cosd(Lambda)] + [e*chord(1) 0];
    geom.pctrl(Nhalf+i2,3) = -gb.pctrl(i2,1).*sind(dtheta) + g0.pctrl(Nhalf+i2,3).*cosd(dtheta) + dw;
    geom.pctrl(1:Nhalf,3) = flipud(geom.pctrl(Nhalf+1:N,3));

    % Evaluate deformed configuration
    [z,t] = LLT(geom,qinf,alpha,e);
    [Pben,Ptor] = fastritz(5,geom.vertex(Nhalf+1:N+1,2)/cosd(Lambda),z,t,Lambda,structprop,method);
    x = [polyval(Pben,b/2/cosd(Lambda));polyval(Ptor,b/2/cosd(Lambda))];
end

cmap = parula(iter);
for i = 1:iter
    ax1.Children(i).Color = cmap(i,:);
    ax2.Children(i).Color = cmap(i,:);
end
clear
close all

set(0,'defaultAxesFontName','Times')
set(0,'defaultAxesFontSize',12)

%% Inputs
method = 'ritz';
Nmodes = 5;

structprop.EI = 380;
structprop.GJ = 143;
e = 0.15;

qinf = 20;
alpha = 2;

Nhalf = 8;                          % number of horseshoe vortices in a semispan
Lambda = 45;                        % c/4 sweep angle (deg)
b0 = 9.01;                          % physical object length
phi = 0;                            % dihedral angle (deg)
ys = linspace(0,1,Nhalf+1).';       % non-dimensional semispan coordinate of legs
twist = linspace(0,0,Nhalf+1).';    % twist distribution on ys (deg)
chord = 5.9/12 + zeros(Nhalf+1,1);  % physical object width distribution on ys

Au = [0.20217 0.17506 0.19269 0.15789 0.16729 0.16283]; % upper surface weights
Al = -Au;                                               % lower surface weights
g0.af = CST_airfoil(Au,Al,51);

%% Process initial geometry
N = 2*Nhalf;
b = b0*cosd(Lambda);        % wingspan
chord = chord/cosd(Lambda); % spanwise chord distribution

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
xlabel('$y/b$','Interpreter','latex','FontSize',14)
ylabel('$\Delta w$ (ft)','Interpreter','latex','FontSize',14)

figure
ax2 = gca;
hold on
xlabel('$y/b$','Interpreter','latex','FontSize',14)
ylabel('$\Delta \theta$ (deg)','Interpreter','latex','FontSize',14)

maxiter = 50; % maximum allowed iterations before terminating

iter = 0;       % iteration counter
geom = g0;      % set iterable to initial state
xLast = [0;0];  % storage for state of previous iteration 

% Perform iteration 0
[z,t] = LLT(geom,qinf,alpha,e);

%% Excursion for extra plots
figure
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
ax3 = nexttile; hold on
ylabel('$\Delta w$ (ft)','Interpreter','latex','FontSize',14)
ax4 = nexttile; hold on
ylabel('$\Delta \theta$ (deg)','Interpreter','latex','FontSize',14)
xlabel(get(gcf,'Children'),'$y/b$','Interpreter','latex','FontSize',14)

for i = [1 2 5 8]
    [Pben,Ptor] = fastritz(i,geom.vertex(Nhalf+1:N+1,2)/cosd(Lambda),z,t,Lambda,structprop,method);
    plot(ax3,ys,polyval(Pben,max(0,gb.vertex(:,2))),'DisplayName',sprintf('%d modes\\quad',i))
    plot(ax4,ys,polyval(Ptor,max(0,gb.vertex(:,2))))
end
ax3.Children(1).DisplayName = regexprep(ax3.Children(1).DisplayName,'\\quad','');
legend(ax3,'Orientation','horizontal','Interpreter','latex','FontSize',12)
ax3.Legend.Layout.Tile = 'north';
ax3.Parent.Parent.Position(3:4) = [900 400];

%% Continue
[Pben,Ptor] = fastritz(Nmodes,geom.vertex(Nhalf+1:N+1,2)/cosd(Lambda),z,t,Lambda,structprop,method);
x = [polyval(Pben,b/2/cosd(Lambda));polyval(Ptor,b/2/cosd(Lambda))];

while (abs(x(1)-xLast(1)) > 1e-3) || (abs(x(2)-xLast(2)) > 1e-2)
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
    [Pben,Ptor] = fastritz(Nmodes,geom.vertex(Nhalf+1:N+1,2)/cosd(Lambda),z,t,Lambda,structprop,method);
    x = [polyval(Pben,b/2/cosd(Lambda));polyval(Ptor,b/2/cosd(Lambda))];
end

cmap = parula(iter);
for i = 1:iter
    ax1.Children(i).Color = cmap(i,:);
    ax2.Children(i).Color = cmap(i,:);
end
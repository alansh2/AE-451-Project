clear
close all

%%
structprop.EI = 8830/144;
structprop.GJ = 13330/144;
e = 0.25;

qinf = 30;
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
geom.af = CST_airfoil(Au,Al,51);

[geom.vertex,geom.pctrl,geom.cctrl] = geom2grid(b,chord,Lambda,phi,twist,ys);

%%
figure
ax1 = gca;
hold on

figure
ax2 = gca;
hold on

N = 2*Nhalf;
cmap = parula(10);
for iter = 1:10
    [z,t] = LLT(geom,qinf,alpha,e);

    [Pben,Ptor] = fastritz(5,geom.vertex(Nhalf+1:N+1,2)/cosd(Lambda),z,t,Lambda,structprop,'ritz');

    % Deform the geometry
    % Regenerate undeformed geometry
    [geom.vertex,geom.pctrl,geom.cctrl] = geom2grid(b,chord,Lambda,phi,twist,ys);
    % Move vertices
    gb = (geom.vertex(Nhalf+1:N+1,1:2)-[e*chord(1) 0])*[cosd(Lambda) sind(Lambda);-sind(Lambda) cosd(Lambda)];
    for i = 1:Nhalf+1
        if gb(i,2) > 0
            dw = polyval(Pben,gb(i,2));
            dtheta = polyval(Ptor,gb(i,2));
            xzp = [gb(i,1) geom.vertex(Nhalf+i,3)]*[cosd(dtheta) -sind(dtheta);sind(dtheta) cosd(dtheta)];
            gb(i,1) = xzp(1);
            geom.vertex(Nhalf+i,3) = xzp(2) + dw;
        end
    end
    geom.vertex(Nhalf+1:N+1,1:2) = gb*[cosd(Lambda) -sind(Lambda);sind(Lambda) cosd(Lambda)]+[e*chord(1) 0];
    geom.vertex(1:Nhalf,3) = flipud(geom.vertex(Nhalf+2:N+1,3));

    % Move control points
    gb = (geom.pctrl(Nhalf+1:N,1:2)-[e*chord(1) 0])*[cosd(Lambda) sind(Lambda);-sind(Lambda) cosd(Lambda)];
    for i = 1:Nhalf
        if gb(i,2) > 0
            dw = polyval(Pben,gb(i,2));
            dtheta = polyval(Ptor,gb(i,2));
            xzp = [gb(i,1) geom.pctrl(Nhalf+i,3)]*[cosd(dtheta) -sind(dtheta);sind(dtheta) cosd(dtheta)];
            gb(i,1) = xzp(1);
            geom.pctrl(Nhalf+i,3) = xzp(2) + dw;
        end
    end
    geom.pctrl(Nhalf+1:N,1:2) = gb*[cosd(Lambda) -sind(Lambda);sind(Lambda) cosd(Lambda)]+[e*chord(1) 0];
    geom.pctrl(1:Nhalf,3) = flipud(geom.pctrl(Nhalf+1:N,3));

    plot(ax1,geom.vertex(Nhalf+1:N+1,2),geom.vertex(Nhalf+1:N+1,3),'Color',cmap(iter,:))
    plot(ax2,geom.pctrl(Nhalf+1:N,2),geom.pctrl(Nhalf+1:N,3),'Color',cmap(iter,:))

    % geom.vertex(Nhalf+1:N+1,3) = geom.vertex(Nhalf+1:N+1,3) + ...
    %     polyval(Pb,max(0,e*chord*sind(Lambda)+geom.vertex(Nhalf+1:N+1,2)/cosd(Lambda)));
    % geom.vertex(1:Nhalf,3) = flipud(geom.vertex(Nhalf+2:N+1,3));
    % geom.pctrl(:,3) = geom.pctrl(:,3) + polyval(Pb,max(0,abs(geom.pctrl(:,2))/cosd(Lambda) - ...
    %     (e*geom.cctrl-(geom.pctrl(:,1)-0.5*(geom.vertex(1:N,1)+geom.vertex(2:N+1,1))))*sind(Lambda)));
end
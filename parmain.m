clear
close all

p = gcp('nocreate');
if isempty(p)
    pool = parpool('Processes');
end

%% Inputs
method = 'lstsq';

structprop.EI = 380;
structprop.GJ = 143;
e = 0.15;

Qinf = [18 23 31.5 51 80];        % freestream dynamic pressure upper limit
alpha = 2;

Nhalf = 8;                      % number of horseshoe vortices in a semispan
b0 = 9.01;                      % wingspan
Lam = [-60 -45 -30 -15 0];      % c/4 sweep angle (deg)
phi = 0;                        % dihedral angle (deg)
ys = linspace(0,1,Nhalf+1).';   % non-dimensional semispan coordinate of legs
twist = zeros(Nhalf+1,1);       % twist distribution on ys (deg)

Au = [0.20217 0.17506 0.19269 0.15789 0.16729 0.16283]; % upper surface weights
Al = -Au;                                               % lower surface weights
g0.af = CST_airfoil(Au,Al,51);

N = 2*Nhalf;

X = cell(1,length(Lam));
for i = 1:length(Lam)
    Lambda = Lam(i);
    b = b0*cosd(Lambda);
    chord = 5.9/12/cosd(Lambda) + zeros(Nhalf+1,1);

    [g0.vertex,g0.pctrl,g0.cctrl] = geom2grid(b,chord,Lambda,phi,twist,ys);
    
    % Align initial geometry with the elastic axis
    croot = chord(1);
    gb.vertex = (g0.vertex(Nhalf+1:N+1,1:2)-[e*croot 0])*[cosd(Lambda) sind(Lambda);-sind(Lambda) cosd(Lambda)];
    gb.pctrl = (g0.pctrl(Nhalf+1:N,1:2)-[e*croot 0])*[cosd(Lambda) sind(Lambda);-sind(Lambda) cosd(Lambda)];
    i1 = find(gb.vertex(:,2)>0,1):Nhalf+1; % query positive positions along the span
    i2 = find(gb.pctrl(:,2)>0,1):Nhalf;
    
    maxiter = 100; % maximum allowed iterations before terminating
    
    qinf = 10:0.5:Qinf(i);
    Xi = zeros(2,length(qinf));
    parfor j = 1:length(qinf)
        iter = 0;       % iteration counter
        geom = g0;      % set iterable to initial state
        xLast = [0;0];  % storage for state of previous iteration 
        
        % Perform iteration 0
        [z,t] = LLT(geom,qinf(j),alpha,e);
        [Pben,Ptor] = fastritz(5,geom.vertex(Nhalf+1:N+1,2)/cosd(Lambda),z,t,Lambda,structprop,method);
        x = [polyval(Pben,b/2/cosd(Lambda));polyval(Ptor,b/2/cosd(Lambda))];
        
        while (abs(x(1)-xLast(1)) > 1e-3) || (abs(x(2)-xLast(2)) > 1e-2)
            if iter >= maxiter
                warning('Exceeded maximum iterations. Solution did not converge within tolerance.')
                break
            end
            % if xLast(1)/(b/2) > 1
            %     break
            % end
            iter = iter + 1;
            xLast = x;
        
            % Move vertices
            dw = polyval(Pben,gb.vertex(i1,2));
            dtheta = polyval(Ptor,gb.vertex(i1,2));
            geom.vertex(Nhalf+i1,1:2) = [gb.vertex(i1,1).*cosd(dtheta)+g0.vertex(Nhalf+i1,3).*sind(dtheta),...
                gb.vertex(i1,2)]*[cosd(Lambda) -sind(Lambda);sind(Lambda) cosd(Lambda)] + [e*croot 0];
            geom.vertex(Nhalf+i1,3) = -gb.vertex(i1,1).*sind(dtheta) + g0.vertex(Nhalf+i1,3).*cosd(dtheta) + dw;
            geom.vertex(1:Nhalf,3) = flipud(geom.vertex(Nhalf+2:N+1,3));
        
            % plot(ax1,ys,[zeros(i1(1)-1,1);dw])
            % plot(ax2,ys,[zeros(i1(1)-1,1);dtheta])
        
            % Move control points
            dw = polyval(Pben,gb.pctrl(i2,2));
            dtheta = polyval(Ptor,gb.pctrl(i2,2));
            geom.pctrl(Nhalf+i2,1:2) = [gb.pctrl(i2,1).*cosd(dtheta)+g0.pctrl(Nhalf+i2,3).*sind(dtheta),...
                gb.pctrl(i2,2)]*[cosd(Lambda) -sind(Lambda);sind(Lambda) cosd(Lambda)] + [e*croot 0];
            geom.pctrl(Nhalf+i2,3) = -gb.pctrl(i2,1).*sind(dtheta) + g0.pctrl(Nhalf+i2,3).*cosd(dtheta) + dw;
            geom.pctrl(1:Nhalf,3) = flipud(geom.pctrl(Nhalf+1:N,3));
        
            % Evaluate deformed configuration
            [z,t] = LLT(geom,qinf(j),alpha,e);
            [Pben,Ptor] = fastritz(5,geom.vertex(Nhalf+1:N+1,2)/cosd(Lambda),z,t,Lambda,structprop,method);
            x = [polyval(Pben,b/2/cosd(Lambda));polyval(Ptor,b/2/cosd(Lambda))];
        end
        Xi(:,j) = x;
    end
    X{i} = Xi;
end

%%
figure
hold on
for i = 1:length(Lam)
    plot(10:0.5:Qinf(i),X{i}(1,:)/(b0*cosd(Lam(i))/2),'DisplayName',sprintf('$\\Lambda=%g^\\circ$',Lam(i)))
end
xlabel('$q_\infty$ (lb/ft$^2$)','Interpreter','latex','FontSize',14)
ylabel('$\Delta w(b)/b$','Interpreter','latex','FontSize',14)
legend('Location','northeast','Interpreter','latex','FontSize',12)
ylim([0 1])
%%
figure
hold on
for i = 1:length(Lam)
    plot(10:0.5:Qinf(i),X{i}(2,:),'DisplayName',sprintf('$\\Lambda=%g^\\circ$',Lam(i)))
end
xlabel('$q_\infty$ (lb/ft$^2$)','Interpreter','latex','FontSize',14)
ylabel('$\Delta \theta(b)$ (deg)','Interpreter','latex','FontSize',14)
legend('Location','northeast','Interpreter','latex','FontSize',12)
ylim([0 5])
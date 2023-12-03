function [Gsol,E] = solvegamma(G0,alpha,vij,un,ua,zeta,af)
% Inputs:
%      G0: initial guess for Gamma/Vinf distribution
%   alpha: angle of attack (deg)
%     vij: LLT influence matrix
%      un: panel normal direction vector
%      ua: panel chordwise direction vector
%    zeta: dl/dA
%      af: airfoil coordinates as a 2-column vector
% Outputs:
%    Gsol: computed solution for Gamma/Vinf
%       E: magnitude of the residual

N = length(G0)/2;
uinf = [cosd(alpha) 0 sind(alpha)];

[Gsol,E] = fminunc(@objfun,G0(1:N));
Gsol = [Gsol;flipud(Gsol)];

    function fval = objfun(Ghalf)
        G = [Ghalf;flipud(Ghalf)];
        v = zeros(N,3);
        for i = 1:3
            v(:,i) = uinf(i) + vij(1:N,:,i)*G/(4*pi);
        end
        alf = atan2d(dot(v,un(1:N,:),2),dot(v,ua(1:N,:),2));
        Cl = zeros(N,1);
        for i = 1:N
            Cl(i) = Panel2D(af,alf(i));
        end
        fval = mean((2*vecnorm(cross(v,zeta(1:N,:),2),2,2).*Ghalf - Cl).^2);
        % fval = max(abs(2*vecnorm(cross(v,zeta(1:N,:),2),2,2).*Ghalf - Cl));
    end
end
% function shape = getshape(twist_ord,bend_ord,bc)

%Inputs:    twist_ord = any positive real integer
%           bend_ord = any positive real integer
%           twist_bc = symbolic array of boundary conditions
%           bend_bc = symbolic array of boundary conditions

if twist_ord ~= height(twist_bc) || bend_ord ~= height(bend_bc)

    fprintf("More Boundary Conditions Required \n")
    return

end

syms y b
syms c [1 twist_ord+1]
syms l [1 bend_ord+1]
w = 0;
twist_coeff = zeros(twist_ord+1);
bend_coeff = zeros(bend_ord+1);
E = 1;
I = 1;
G = 1;
J = 1;

theta = poly2sym(flip(c),y);

% for i = 1:(twist_ord+1) 
% 
%     theta = theta + twist_coeff(i)*y^(i-1);
% 
% end

if max(twist_bc(:,3)) > 0
    for j = 1:max(twist_bc(:,3))
        theta(j+1) = diff(theta,y,j);
    end
else
    fprintf("Error")
    return
end


for i = 1:height(twist_bc)

    twist_cond(i) = theta(twist_bc(i,3)+1) ==

end

% thetap = diff(theta,y);
% c0 = solve(subs(theta,y,0) == 0, twist_coeff(1));
% c1 = solve(subs(thetap,y,b) == 0, twist_coeff(2));
% theta = subs(theta,twist_coeff(1),c0);
% theta = subs(theta,twist_coeff(2),c1);
% thetap = diff(theta,y);
% 
% 
% 
% 
% for i = 1:(bend_ord+1) 
% 
%     bend_coeff(i) = str2sym("l" + num2str(i-1));
%     w = w + bend_coeff(i)*y^(i-1);
% 
% end
% 
% 
% wp = diff(w,y);
% wpp = diff(wp,y);
% l0 = solve(subs(w,y,0) == 0, bend_coeff(1));
% l1 = solve(subs(wp,y,0) == 0, bend_coeff(2));
% l2 = solve(subs(wpp,y,b),bend_coeff(3));
% w = subs(w,bend_coeff(1),l0);
% w = subs(w,bend_coeff(2),l1);
% w = subs(w,bend_coeff(3),l2);
% wp = diff(w,y);
% wpp = diff(wp,y);
% 
% U = (1/2)*int(E*I*wpp^2 + G*J*thetap^2,0,b);


% N = size(dl,1);
% dybar = vecnorm(dl(:,1:2),2,2);
% ybar = [0;cumsum(dybar(1:N-1))] + dybar/2;
% bbar = ybar(N) + dybar(N)/2;
% l3 = -sum(Lbar.*(ybar.^3-3*bbar*ybar.^2).*dybar)/(12*E*I*bbar^2);
% l2 = -3*l3*bbar;
% polyout = [l3 l2 0 0];
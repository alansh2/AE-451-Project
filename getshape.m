% function shape = getshape(twist_ord,bend_ord,twist_bc,bend_bc)

%Inputs:    twist_ord = any positive real integer
%           bend_ord = any positive real integer
%           twist_bc = symbolic array of boundary conditions [set to,
%           y-value, differential order]
%           bend_bc = symbolic array of boundary conditions [set to,
%           y-value, differential order]

syms y b                            %Move down to variable initialization when this becomes a function
twist_ord = 2;                      %Using to run as a stand alone script
bend_ord = 3;                       %Using to run as a stand alone script
twist_bc = [0 0 0; 0 b 1];          %Using to run as a stand alone script
bend_bc = [0 0 0; 0 0 1; 0 b 2];    %Using to run as a stand alone script

%% Boundary Condition Check

if twist_ord ~= height(twist_bc) || bend_ord ~= height(bend_bc) %Make sure there are enough B.C.s to solve the system

    fprintf("More Boundary Conditions Required \n") %If not enough B.C.s, output error
    return                                          %and end function.

end

%% Variable Initialization

syms c [1 twist_ord+1]  %Initialize array of syms for twist shape function
c = flip(c);            %Reverse the order of variables such that poly2sim puts the smallest "name" next to the lowest order variable
syms l [1 bend_ord+1]   %Initialize array of syms for bending shape function
l = flip(l);            %Reverse the order of variables such that poly2sim puts the smallest "name" next to the lowest order variable

%Structural Properties
E = 1;
I = 1;
G = 1;
J = 1;

theta = poly2sym(c,y);  %Create polynomial of order equal to length of c with coefficients c and variable y
w = poly2sym(l,y);      %Create polynomial of order equal to length of l with coefficients l and variable y

%% Shape Function Differentiation

if max(twist_bc(:,3)) > 0                   %Find the highest order derivative necessary for twist
    for j = 1:max(twist_bc(:,3))            %Iterate to get each required derivative
        theta(j+1) = diff(theta(1),y,j);
    end
else
    fprintf("Error")                        %Output error if there are no derivative boundary conditions
    return
end

if max(bend_bc(:,3)) > 0                    %Find the highest order derivative necessary for bending
    for j = 1:max(bend_bc(:,3))             %Iterate to get each required derivative
        w(j+1) = diff(w(1),y,j);
    end
else
    fprintf("Error")                        %Output error if there are no derivative boundary conditions
    return
end

%% System of Equations Solver

%Initialize variables to store new expressions for each coefficient
c_new = c;
l_new = l;

%Solution for twist
for i = 1:height(twist_bc)

    twist_cond(i) = subs(theta(twist_bc(i,3)+1),y,twist_bc(i,2)) == twist_bc(i,1);  %Determine the B.C. equation when values have been substituted

    for j = 1:length(c)
        check = has(twist_cond(i),c(j));                    %See what the highest number name coefficient contained in the i^th expression is
        if check == true
            if j == length(c)
                c_new(j) = solve(twist_cond(i),c(j));       %If only the last ("constant") coefficient is contained, solve for it
            elseif has(twist_cond(i),c(j+1)) == true
                c_new(j+1) = solve(twist_cond(i),c(j+1));   %If there is an expression with more than one coefficient, solve for the lower of them
            else
                c_new(j) = solve(twist_cond(i),c(j));       %If there is only a single coefficient, solve for it directly
            end
            break                                           %End this boundary condition
        else
            continue                                        %Check for the next coefficient
        end
    end

end

for i = 1:length(c)                             %Iterate to substitute in the coefficients such that the equation is in terms of just one
    theta(1) = subs(theta(1),c(i),c_new(i));
end

%Solution for bending
for i = 1:height(bend_bc)

    bend_cond(i) = subs(w(bend_bc(i,3)+1),y,bend_bc(i,2)) == bend_bc(i,1);

    for j = 1:length(l)
        check = has(bend_cond(i),l(j));                     %See what the highest number name coefficient contained in the i^th expression is
        if check == true
            if j == length(l)
                l_new(j) = solve(bend_cond(i),l(j));        %If only the last ("constant") coefficient is contained, solve for it
            elseif has(bend_cond(i),l(j+1)) == true
                l_new(j+1) = solve(bend_cond(i),l(j+1));    %If there is an expression with more than one coefficient, solve for the lower of them
            else
                l_new(j) = solve(bend_cond(i),l(j));        %If there is only a single coefficient, solve for it directly
            end
            break                                           %End this boundary condition
        else
            continue                                        %Check for the next coefficient
        end
    end

end

for i = 1:length(l)                     %Iterate to substitute in the coefficients such that the equation is in terms of just one
    w(1) = subs(w(1),l(i),l_new(i));
end

%% Shape Function Differentiation 2

theta(2) = diff(theta(1),y);    %Find the 1st derivative of twist shape
w(2) = diff(w(1),y);            %Find the 1st derivative of bending shape
w(3) = diff(w(2),y);            %Find the 2nd derivative of bending shape

%% Potential Energy

U = (1/2)*int(E*I*w(3)^2 + G*J*theta(2)^2,0,b); %Integrate using the potential energy definition


% N = size(dl,1);
% dybar = vecnorm(dl(:,1:2),2,2);
% ybar = [0;cumsum(dybar(1:N-1))] + dybar/2;
% bbar = ybar(N) + dybar(N)/2;
% l3 = -sum(Lbar.*(ybar.^3-3*bbar*ybar.^2).*dybar)/(12*E*I*bbar^2);
% l2 = -3*l3*bbar;
% polyout = [l3 l2 0 0];












%% Old Code


% line 37
% for i = 1:(twist_ord+1) 
% 
%     theta = theta + twist_coeff(i)*y^(i-1);
% 
% end

% for i = 1:length(twist_cond)
% 
%     c(i) = solve()
% 
% end


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
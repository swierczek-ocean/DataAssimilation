function [X_t,M] = L96_R_function(X,ode_rhs_fun,gradient_fun,dt,jump,n)
%% R function calculation for 4DVar on system
% Calculates trajectory of system for jump steps.
% The first four are performed with RK4, then the remainder
% are done with AB4. So jump must be greater than 4.
% Also calculates derivate of this algorithm.
%%

%% preliminaries
FEvals = [ode_rhs_fun(X),zeros(n,3)];
dt = dt/4;
M = eye(n);
dX = zeros(n,n,jump);
dX(:,:,1) = M;
GFEvals = zeros(n,n,jump);
GFEvals(:,:,1) = gradient_fun(X);
%%

for ii=2:4
    %% RK4
    k1=ode_rhs_fun(X);
    k2=ode_rhs_fun(X+0.5*dt.*k1);
    k3=ode_rhs_fun(X+0.5*dt.*k2);
    k4=ode_rhs_fun(X+dt.*k3);       
    %%
    
    %% derivative of RK4
    p1 = gradient_fun(X+0.5*dt.*k1);
    p2 = p1.*(eye(n)+0.5*dt.*gradient_fun(X));
    p3 = gradient_fun(X+0.5*dt.*k2);
    p4 = p3.*(eye(n)+0.5*dt.*p2);
    p5 = gradient_fun(X+0.5.*k3);
    p6 = p5.*(eye(n)+0.5*dt.*p4);       
    
    dX(:,:,ii) = (eye(n) + (1/6)*dt.*gradient_fun(X) + ...
        (1/3)*dt.*(p2 + p4) + (1/6)*dt.*p6)*dX(:,:,ii-1);
    %%
    
    X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);
    FEvals(:,ii) = ode_rhs_fun(X);
    GFEvals(:,:,ii) = gradient_fun(X);
    
end

for ii=5:jump
    %% AB4
    Temp_1 = (55/24).*FEvals(:,4) - ...
        (59/24).*FEvals(:,3) + ...
        (37/24).*FEvals(:,2) - ...
        (3/8).*FEvals(:,1);
    X = X + dt.*Temp_1;
    Temp_2 = ode_rhs_fun(X);
    FEvals(:,1:3) = FEvals(:,2:4);
    FEvals(:,4) = Temp_2;
    %%
    
    %% derivative of AB4
    dX(:,:,ii) = dX(:,:,ii-1) + ...
        dt.*((55/24).*GFEvals(:,:,ii-1)*dX(:,:,ii-1) - ...
        (59/24).*GFEvals(:,:,ii-2)*dX(:,:,ii-2) + ...
        (37/24).*GFEvals(:,:,ii-3)*dX(:,:,ii-3) - ...
        (3/8).*GFEvals(:,:,ii-4)*dX(:,:,ii-4));
    %%
    
    GFEvals(:,:,ii) = gradient_fun(X);
end

X_t = X;
M = dX(:,:,end);
end


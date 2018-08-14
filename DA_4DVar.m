function [X_star_t,X_star,Time_Series,J,Cov] = DA_4DVar(X,ode_rhs_fun,gradient_fun,Cov,H,mu,dt,jump,Obs,ObsVar,n)
%% 4DVar
% Performs one step of a 4DVar DA scheme
% First, a function R is specified so that we can perform
% a least squares minimization on the function R'*R
%%

fun = @(x)ACC_4DVar_R_wrapper(x,ode_rhs_fun,gradient_fun,Cov,H,mu,dt,jump,Obs,ObsVar,n);

options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
    'SpecifyObjectiveGradient',true,...
    'Diagnostics','off',...
    'Display','off',...         % 'iter-detailed'
    'CheckGradient',false);

[X_star,~,~,~,~,~,J] = lsqnonlin(fun,X,[],[],options);

Time_Series = [X_star,zeros(n,jump-1)];

%% preliminaries
FEvals = [ode_rhs_fun(X_star),zeros(n,3)];
%%

X = X_star;

for ii=2:4
    %% RK4
    k1=ode_rhs_fun(X);
    k2=ode_rhs_fun(X+0.5*dt.*k1);
    k3=ode_rhs_fun(X+0.5*dt.*k2);
    k4=ode_rhs_fun(X+dt.*k3);       
    X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);
    FEvals(:,ii) = ode_rhs_fun(X);
    Time_Series(:,ii) = X;
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
    Time_Series(:,ii) = X;
    %%
end

X_star_t = X;

[~,M] = ACC_R_function(X_star,ode_rhs_fun,gradient_fun,dt,jump,n);
% Cov = M\(2*(J'*J))*(M');
end


function [X,FEvals] = ODE_AB4_auto(X,FEvals,ode_rhs_fun,dt)
%% Fourth Order Adams-Bashforth Method
% Input FEval must be a matrix with 4 columns.
% Each column represents the function evaluated at a specific time.
% The rightmost column is at the corrent time, the 2nd from right
% is the previous time step, etc.
% Described in https://en.wikipedia.org/wiki/Linear_multistep_method
% This version is for autonomous ODEs dx/dt = F(x)
%%

%% AB4
Temp_1 = (55/24).*FEvals(:,4) - ...
    (59/24).*FEvals(:,3) + ...
    (37/24).*FEvals(:,2) - ...
    (3/8).*FEvals(:,1);
X = X + dt.*Temp_1;
%%

%% Function Evaluation
Temp_2 = ode_rhs_fun(X);
%%
%% Configure Output
FEvals(:,1:3) = FEvals(:,2:4);
FEvals(:,4) = Temp_2;
%%

end


function [X,FEvals] = ODE_ABMPC4(X,FEvals,ode_rhs_fun,dt,time_current)
%% Fourth Order Adams-Bashforth / Adams-Moulton Predictor-Corrector Method
% Input X must be a matrix with 4 columns.
% Each column represents the state vector at a specific time.
% The rightmost column is the corrent time, the 2nd from right
% is the previous time step, etc.
% Described in https://en.wikipedia.org/wiki/Linear_multistep_method
% This version is for time dependent ODEs dx/dt = F(t,x)
%%

%% Predictor
Temp_1 = (55/24).*FEvals(:,4) - ...
    (59/24).*FEvals(:,3) + ...
    (37/24).*FEvals(:,2) - ...
    (3/8).*FEvals(:,1);
Pred = X + dt.*Temp_1;
%%

%% Corrector
Temp_2 = ode_rhs_fun(time_current+dt,Pred);
Temp_3 = (251/720).*Temp_2 + ...
    (646/720).*FEvals(:,4) - ...
    (264/720).*FEvals(:,3) + ...
    (106/720).*FEvals(:,2) - ...
    (19/720).*FEvals(:,1);
X = X + dt.*Temp_3;
%%
%% Configure Output
FEvals(:,1:3) = FEvals(:,2:4);
FEvals(:,4) = Temp_2;
%%

end


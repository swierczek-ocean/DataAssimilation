function [X,FEvals] = ODE_RK4_auto_start(ode_rhs_fun,init_cond,dt)
%% Four step seed for P-C method using RK4
% This uses a 4th order Runge-Kutta ODE solver to seed 
% a fourth order Adams-Bashforth / Adams-Moulton predictor-corrector
% where the ODE has an autonomous RHS
%%
n = size(init_cond,1);      %% the state vector needs to be a column
X = init_cond;
FEvals = [ode_rhs_fun(init_cond),zeros(n,3)];

for ii=2:4
   X = ODE_RK4_auto(X,ode_rhs_fun,dt); 
   FEvals(:,ii) = ode_rhs_fun(X);
end
end


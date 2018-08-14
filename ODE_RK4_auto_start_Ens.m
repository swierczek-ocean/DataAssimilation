function [Ensemble,EnFEvals] = ODE_RK4_auto_start_Ens(ode_rhs_fun,Ensemble,dt)
%% Four step seed for P-C method using RK4
% This uses a 4th order Runge-Kutta ODE solver to seed 
% a fourth order Adams-Bashforth / Adams-Moulton predictor-corrector
% where the ODE has an autonomous RHS
%%
[n,Ne] = size(Ensemble);
EnFEvals = zeros(n,4,Ne);
EnFEvals(:,1,:) = ode_rhs_fun(Ensemble);

for ii=2:4
    Ensemble = ODE_RK4_auto(Ensemble,ode_rhs_fun,dt);
    EnFEvals(:,ii,:) = ode_rhs_fun(Ensemble);
end

end


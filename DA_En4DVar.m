function [Ensemble,X_star_t,spread,Time_Series,Cov] = DA_En4DVar(X_star_t,Ensemble,ode_rhs_fun,gradient_fun,H,Cov,mu,dt,jump,Obs,ObsVar,n,L,alpha)
%% Ensemble Kalman Filter with 4DVar
% Hybrid DA algorithm:
% 1) Takes mean of ensemble at current time, performs 4DVar.
% 2) Performs EnKF on rest of ensemble
% 3) Recenters output ensemble about X_star_t
%%

%% 4DVar
[X_star_t,~,Time_Series,~,~] = DA_4DVar(X_star_t,ode_rhs_fun,...
    gradient_fun,Cov,H,mu,dt,jump,Obs,ObsVar,n);
%%

%% EnKF
[Ensemble,mu_a,spread,Cov] = DA_SqEnKF_plus(Ensemble,ode_rhs_fun,H,jump,Obs,ObsVar,L,alpha,dt);
%%

%% Recentering
Ensemble = Ensemble - mu_a + X_star_t;
%%

end


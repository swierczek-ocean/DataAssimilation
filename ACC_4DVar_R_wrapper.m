function [R,Jacobian] = ACC_4DVar_R_wrapper(X,ode_rhs_fun,gradient_fun,Cov,H,mu,dt,jump,Obs,ObsVar,n)
%% Objective Function for 4DVar
% Creates function to be minimized for 4DVar DA
%%

%% Run model
% call outside function that generates time series, gradients,
% and derivatives.
%%
[X_t,M] = ACC_R_function(X,ode_rhs_fun,gradient_fun,dt,jump,n);
%%

%% Create R vector and Jacobian
R = [real(sqrtm(Cov))\(X-mu);(H*X_t-Obs)./sqrt(ObsVar)];
Jacobian = [real(sqrtm(Cov))\eye(n);H*M./sqrt(ObsVar)];
%%

end


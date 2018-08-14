function [X_star_t,spread,Cov,Time_Series] = DA_EDA(X,ode_rhs_fun,gradient_fun,Cov,H,mu,dt,jump,Obs,ObsVar,n,L,alpha,Ne)
%% Ensemble 4DVar / RTO algorithm
% Ensemble variational algorithm:
% 1) Takes mean of ensemble at current time, performs 4DVar.
% 2) Performs 4DVar on perturbations of ensemble mean
%%

%% Preliminaries
m = size(H,1);
Cov = (1+alpha).*L.*Cov;
Cov = 0.5*(Cov+Cov');
mu_perts = sqrtm(Cov)*randn(n,Ne-1);
Obs_perts = sqrt(ObsVar)*randn(m,Ne-1);

%%

%% 4DVar + EDA

[X_star_t,~,Time_Series,~,~] = DA_4DVar(X,ode_rhs_fun,gradient_fun,...
    Cov,H,mu,dt,jump,Obs,ObsVar,n);         % unperturbed 4DVar
Ensemble = [X_star_t,zeros(n,Ne-1)];

for ii=1:Ne-1
    Ensemble(:,ii+1) = DA_4DVar(X,ode_rhs_fun,gradient_fun,...
        Cov,H,mu+mu_perts(:,ii),dt,jump,Obs+Obs_perts(:,ii),ObsVar,n);
end
%%

%% Recentering
Ensemble(:,1:Ne) = Ensemble(:,1:Ne) - mean(Ensemble(1:Ne),2) + X_star_t;
%%

%% Additional Outputs
Cov = cov(Ensemble');
spread = sqrt(trace(Cov)/n);
%%

end


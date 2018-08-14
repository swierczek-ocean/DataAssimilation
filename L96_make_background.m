function [mu,BCov] = L96_make_background(fun,dt,long_time,n)
%% Makes an initial ensemble for L96 DA
% runs the L96 from a random initial condition for a 
% long time, and grabs Ne (semi)random states (see implementation). 
% These become our starting ensemble.
%%

%% preliminaries
IC = unifrnd(-1,1,n,1);
num_steps = floor(long_time/dt);
init_steps = floor(0.2*num_steps);
keep_steps = num_steps-init_steps;
Ensemble = zeros(n,keep_steps);
%%

%% make ensemble
[X,FEvals] = ODE_RK4_auto_start(fun,IC,dt);


for ii=1:init_steps
    [X,FEvals] = ODE_AB4_auto(X,FEvals,fun,dt);  
end

Ensemble(:,1) = X;

for jj=1:keep_steps-1
    [Ensemble(:,jj+1),FEvals] = ODE_AB4_auto(Ensemble(:,jj),FEvals,fun,dt);
end

mu = mean(Ensemble,2);
BCov = cov(Ensemble');

end


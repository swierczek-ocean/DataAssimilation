function Ensemble = L96_make_ensemble(fun,Ne,dt,long_time,n)
%% Makes an initial ensemble for L96 DA
% runs the L96 from a random initial condition for a 
% long time, and grabs Ne (semi)random states (see implementation). 
% These become our starting ensemble.
%%

%% preliminaries
IC = unifrnd(-1,1,n,1);
num_steps = floor(long_time/dt);
Ensemble = zeros(n,Ne);
indices = randperm(num_steps-1000,Ne)+1000; % don't want early iterations
indices = sort(indices);
%%

%% make ensemble
[X,FEvals] = ODE_RK4_auto_start(fun,IC,dt);
kk = 1;

for ii=1:Ne
   while kk<indices(ii)
      [X,FEvals] = ODE_AB4_auto(X,FEvals,fun,dt); 
      kk = kk + 1;
   end
   Ensemble(:,ii) = X;
end
end


tic()           % timer
clc
clear
close all

%% preliminaries
ACC_Colors
n = 40;             % dimension of L96 system
sqn = sqrt(n);
Ne = 20;            % ensemble size
spinup_time = 100;  % for getting onto attractor
exp_time = 1;       % dimensionless time units of DA experiment
long_time = 1000;   % long simulation for creating initial ensemble
dt = 0.01;          % model time step
jump = 10;          % number of model time steps between observations
k = 2;              % observe every kth state variable
F = 8*ones(n,1);    % free parameter on L96 RHS (F = 8 leads to chaotic solutions)
r = 5;            % localization radius
alpha = 0.1;       % ensemble inflation parameter
ObsVar = 1;         % measurement/observation variance
beta = 0.5;
color1 = 11;
color2 = 21;
spinup_iter = floor(spinup_time/dt);    % number of spinup model time steps
exp_iter = floor(exp_time/dt);          % number of experiment model time steps
q = floor(exp_iter/jump);               % number of observed time steps
ObsTimes = jump+1:jump:(exp_iter+jump); % vector of times when observation occurs
%%

%% setup & utilities
[L1,L2] = L96_get_matrices(n);          % makes matrices for matrix-vector execution of L96
[H,m] = L96_get_H(n,k);                 % creates observation operator
L96fun = @(x)((L1*x).*(L2*x) - x + F);  % Lorenz '96 dynamical system
gradient_fun = @(x)L96_gradient(x,L1,L2,n);     % Lorenz '96 gradient
x_start = unifrnd(-1,1,n,1);            % random initial condition
L = ACC_Localize(n,r);                  % localization matrix for covariance
%%

%% spinup for initial conditions
% Run a long simulation to get from initial condition 
% onto L96 attractor. Don't save anything except final time step.
% We use a fourth order Adams-Bashforth linear multistep method.
% This requires a fourth order Runge-Kutta method to get started.
% (The 'auto' in the function names refers to the L96 ODE being autonomous)
%%
[X,FEvals] = ODE_RK4_auto_start(L96fun,x_start,dt);

for ii=1:spinup_iter
    [X,FEvals] = ODE_AB4_auto(X,FEvals,L96fun,dt);
end
%%

%% Make background values
[mu,BCov] = L96_make_background(L96fun,dt,long_time,n);
spread = sqrt(trace(BCov)/n);
init_cond = mu;
%%

%% Make ensemble
Ensemble = L96_make_ensemble(L96fun,Ne,dt,long_time,n);
%%

%% experiment + observations + DA
% Run for exp_time, and observe every kth variable at every jump_th 
% model time step. Perform DA with each set of observations.
%%

[Ensemble,EnFEval] = ODE_RK4_auto_start_Ens(L96fun,Ensemble,dt);
X_star_t_EDA = mean(Ensemble,2);
ErrorVecEDA = zeros(exp_iter,1);
spreadVecEDA = spread.*ones(exp_iter,1);
Time_Series_True = [X,zeros(n,exp_iter)];    % array for storing full true state    
TimeSeriesEDA = zeros(n,exp_iter+1);         % array for storing full DA
total_steps = 0;


num_steps = ObsTimes(1)
for ii=2:num_steps
    [Time_Series_True(:,ii),FEvals] = ODE_AB4_auto(Time_Series_True(:,ii-1),FEvals,L96fun,dt);
end

Obs = H*Time_Series_True(:,num_steps);
[X_star_t_EDA,spread,Cov,Time_Series] = DA_EDA(X_star_t_EDA,L96fun,...
    gradient_fun,BCov,H,mu,dt,jump,Obs,ObsVar,n,L,alpha,Ne);
TimeSeriesEDA(:,1:(num_steps-1)) = Time_Series(:,1:(num_steps-1));
total_steps = total_steps + num_steps;

for kk=2:q
    num_steps = ObsTimes(kk)-ObsTimes(kk-1)
    
    for ii=(ObsTimes(kk-1)+1):ObsTimes(kk)
        [Time_Series_True(:,ii),FEvals] = ODE_AB4_auto(Time_Series_True(:,ii-1),FEvals,L96fun,dt);
    end
    
    Obs = H*Time_Series_True(:,ObsTimes(kk));
    [X_star_t_EDA,spread,Cov,Time_Series] = DA_EDA(X_star_t_EDA,L96fun,...
        gradient_fun,Cov,H,X_star_t_EDA,dt,jump,Obs,ObsVar,n,L,alpha,Ne);
    TimeSeriesEDA(:,ObsTimes(kk-1):(ObsTimes(kk)-1)) = Time_Series(:,1:(num_steps));
    spreadVecEDA(ObsTimes(kk-1):(ObsTimes(kk)-1)) = spread.*ones(num_steps,1);
end

TimeSeriesEDA(:,end) = X_star_t_EDA;

ErrorEDA = TimeSeriesEDA - Time_Series_True;

for ll=1:exp_iter+1
   ErrorVecEDA(ll) = norm(ErrorEDA(:,ll),2)/sqn; 
end

error_parameter = mean(ErrorVecEDA(10*jump:end));
fprintf('Average RMSE: %g\n',error_parameter)
fprintf('Average spread: %g\n',mean(spreadVecEDA(10:end)))

%% error plot
set(gcf, 'Position', [25, 25, 1600, 900])
h1 = plot(spreadVecEDA,'Color',Color(:,color1),'LineWidth',2.2);
hold on
h2 = plot(ErrorVecEDA,'Color',Color(:,color2),'LineWidth',2.2);
title('RMSE & spread for 4DVar')
% axis([0 q 0 3*error_parameter])
xlabel('time')
legend([h1(1),h2(1)],'spread','RMSE')
print('Test_L96_EDA','-djpeg')
hold off
%%

toc()
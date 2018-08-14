tic()           % timer
clc
clear
close all

%% preliminaries
ACC_Colors
n = 40;             % dimension of L96 system
Ne_Sq = 40;         % ensemble size
spinup_time = 100;  % for getting onto attractor
exp_time = 6;       % dimensionless time units of DA experiment
long_time = 1000;   % long simulation for creating initial ensemble
dt = 0.01;          % model time step
jump = 10;          % number of model time steps between observations
k = 2;              % observe every kth state variable
F = 8*ones(n,1);    % free parameter on L96 RHS (F = 8 leads to chaotic solutions)
r1 = 5;             % SqEnKF localization radius
r2 = 5;             % 4DVar localization radius
alpha1 = 0.10;      % SqEnKF inflation parameter
alpha2 = 0.10;      % 4DVar inflation parameter
ObsVar = 1;         % measurement/observation variance
beta = 0.3;
color1 = 21;
color2 = 11;
spinup_iter = floor(spinup_time/dt);    % number of spinup model time steps
exp_iter = floor(exp_time/dt);          % number of experiment model time steps
q = floor(exp_iter/jump);               % number of observed time steps
q_split = ceil(q/2);                    % run EnKF until halfway, then do 4DVar
ObsTimes = jump:jump:(exp_iter+jump); % vector of times when observation occurs
%%

%% setup & utilities
[L1,L2] = L96_get_matrices(n);          % makes matrices for matrix-vector execution of L96
[H,m] = L96_get_H(n,k);                 % creates observation operator
L96fun = @(x)((L1*x).*(L2*x) - x + F);  % Lorenz '96 dynamical system
gradient_fun = @(x)L96_gradient(x,L1,L2,n);     % Lorenz '96 gradient
x_start = unifrnd(-1,1,n,1);            % random initial condition
L_SqEnKF = ACC_Localize(n,r1);          % SqEnKF localization matrix for covariance
L_En4DVar = ACC_Localize(n,r2);         % En4DVar localization matrix for covariance
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

%% Make ensemble
EnsembleSqEnKF = L96_make_ensemble(L96fun,Ne_Sq,dt,long_time,n);
spread = sqrt(trace(cov(EnsembleSqEnKF'))/n);
%%

%% experiment + observations + DA
% Run for exp_time, and observe every kth variable at every jump_th 
% model time step. Perform DA with each set of observations.
%%

[EnsembleSqEnKF,EnFEvalSqEnKF] = ODE_RK4_auto_start_Ens(L96fun,EnsembleSqEnKF,dt);

TimeSeriesEn4DVar = zeros(n,exp_iter);         % array for storing full 4DVar

spreadVecEn4DVar = spread.*ones(1,exp_iter);   

ErrorVecEn4DVar = zeros(1,exp_iter);

Time_Series_True = [X,zeros(n,exp_iter-1)];    % array for storing full true state  
total_steps = 0;

%% from start to first observartions

num_steps = ObsTimes(1)
for ii=2:num_steps
    [Time_Series_True(:,ii),FEvals] = ODE_AB4_auto(Time_Series_True(:,ii-1),FEvals,L96fun,dt);
end

Obs = H*Time_Series_True(:,num_steps);

%% SqEnKF
for jj=2:num_steps
    for mm=1:Ne_Sq
        [EnsembleSqEnKF(:,mm),EnFEvalSqEnKF(:,:,mm)] = ...
            ODE_AB4_auto(EnsembleSqEnKF(:,mm),EnFEvalSqEnKF(:,:,mm),L96fun,dt);
    end
    TimeSeriesEn4DVar(:,jj) = mean(EnsembleSqEnKF,2);
end
[EnsembleSqEnKF,mu_a,spread] = DA_SqEnKF(EnsembleSqEnKF,H,Obs,ObsVar,L_SqEnKF,alpha1);
TimeSeriesEn4DVar(:,num_steps) = mu_a;

%%

total_steps = total_steps + num_steps;
%% 

%% loop for EnKF spinup

for kk=2:q_split  
    num_steps = ObsTimes(kk)-ObsTimes(kk-1)
    
    for ii=(ObsTimes(kk-1)+1):ObsTimes(kk)
        [Time_Series_True(:,ii),FEvals] = ODE_AB4_auto(Time_Series_True(:,ii-1),FEvals,L96fun,dt);
    end
    
    Obs = H*Time_Series_True(:,ObsTimes(kk));
    
    %% SqEnKF
    for jj=1:4
        EnsembleSqEnKF = ODE_RK4_auto(EnsembleSqEnKF,L96fun,dt);
        EnFEvalSqEnKF(:,jj,:) = L96fun(EnsembleSqEnKF);
        TimeSeriesEn4DVar(:,ObsTimes(kk-1)+jj) = mean(EnsembleSqEnKF,2);
    end
    
    for jj=5:num_steps
        for mm=1:Ne_Sq
            [EnsembleSqEnKF(:,mm),EnFEvalSqEnKF(:,:,mm)] = ...
                ODE_AB4_auto(EnsembleSqEnKF(:,mm),EnFEvalSqEnKF(:,:,mm),L96fun,dt);
        end
        TimeSeriesEn4DVar(:,ObsTimes(kk-1)+jj) = mean(EnsembleSqEnKF,2);
    end
    [EnsembleSqEnKF,mu_a,spread,P_a] = DA_SqEnKF(EnsembleSqEnKF,H,Obs,ObsVar,L_SqEnKF,alpha1);
    TimeSeriesEn4DVar(:,ObsTimes(kk)) = mu_a;
    spreadVecEn4DVar(ObsTimes(kk-1):(ObsTimes(kk)-1)) = spread.*ones(num_steps,1);
    %%
end

X_star_t_En4DVar = mu_a;
CovEn4DVar = (1+alpha1)*L_SqEnKF.*P_a;
EnsembleEn4DVar = EnsembleSqEnKF;

for kk=q_split+1:q  
    num_steps = ObsTimes(kk)-ObsTimes(kk-1)
    
    for ii=(ObsTimes(kk-1)+1):ObsTimes(kk)
        [Time_Series_True(:,ii),FEvals] = ODE_AB4_auto(Time_Series_True(:,ii-1),FEvals,L96fun,dt);
    end
    
    Obs = H*Time_Series_True(:,ObsTimes(kk));
    
    %% En4DVar
    [EnsembleEn4DVar,X_star_t_En4DVar,spread,Time_Series,CovEn4DVar] = DA_En4DVar(X_star_t_En4DVar,...
        EnsembleEn4DVar,L96fun,gradient_fun,H,CovEn4DVar,X_star_t_En4DVar,dt,jump,Obs,ObsVar,n,L_En4DVar,alpha1);
    TimeSeriesEn4DVar(:,ObsTimes(kk-1):(ObsTimes(kk)-1)) = Time_Series(:,1:(num_steps));
    spreadVecEn4DVar(ObsTimes(kk-1):(ObsTimes(kk)-1)) = spread.*ones(num_steps,1);
    %%
end
%%

TimeSeriesEn4DVar(:,end) = X_star_t_En4DVar;
spreadVecEn4DVar(:,end) = spread;
ErrorEn4DVar = TimeSeriesEn4DVar - Time_Series_True;

for ll=1:exp_iter
    ErrorVecEn4DVar(ll) = norm(ErrorEn4DVar(:,ll),2); 
end

error_parameter_1 = mean(ErrorVecEn4DVar(10*jump:end));
fprintf('Average RMSE for En4DVar: %g\n',error_parameter_1)
fprintf('Average spread: %g\n',mean(spreadVecEn4DVar(10:end)))

%% error plot
set(gcf, 'Position', [25, 25, 1600, 900])
h1 = plot(ErrorVecEn4DVar,'Color',Color(:,color1),'LineWidth',2.2);
hold on
h2 = plot(spreadVecEn4DVar,'Color',Color(:,color2),'LineWidth',2.2);
title('RMSE & spread')
xlabel('time')
ylabel('RMSE')
legend([h1(1),h2(1)],'4DVar','spread')
print('Test_L96_En4DVar_2','-djpeg')
hold off
%%

toc()
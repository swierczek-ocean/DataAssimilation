tic()           % timer
clc
clear
close all

%% preliminaries
ACC_Colors
n = 40;             % dimension of L96 system
sqn = sqrt(n);
Ne_Sq = 40;         % ensemble size
Ne_En = 30;         % ensemble size
Ne_EDA = 20;        % ensemble size
spinup_time = 100;  % for getting onto attractor
exp_time = 12;      % dimensionless time units of DA experiment
long_time = 1000;   % long simulation for creating initial ensemble
dt = 0.01;          % model time step
jump = 10;          % number of model time steps between observations
k = 2;              % observe every kth state variable
F = 8*ones(n,1);    % free parameter on L96 RHS (F = 8 leads to chaotic solutions)
r1 = 5;             % En4DVar localization radius
r2 = 5;             % 4DVarlocalization radius
r3 = 5.4;           % SqEnKF localization radius
r4 = 4.2;           % EDA localization radius
alpha1 = 0.10;      % En4DVar inflation parameter
alpha2 = 0.10;      % 4DVar inflation parameter
alpha3 = 0.08;      % SqEnKF inflation parameter
alpha4 = 0.05;      % EDA inflation parameter
ObsVar = 1;         % measurement/observation variance
sigma = sqrt(ObsVar);
beta = 0.1;
color1 = 15;
color2 = 11;
color3 = 19;
color4 = 8;
color5 = 9;
spinup_iter = floor(spinup_time/dt);    % number of spinup model time steps
exp_iter = floor(exp_time/dt);          % number of experiment model time steps
q = floor(exp_iter/jump);               % number of observed time steps
q_split = floor((5/6)*q);
ObsTimes = jump+1:jump:(exp_iter+jump); % vector of times when observation occurs
%%

%% setup & utilities
[L1,L2] = L96_get_matrices(n);          % makes matrices for matrix-vector execution of L96
[H,m] = L96_get_H(n,k);                 % creates observation operator
mdim = size(H,1);                       % number of observed state variables
L96fun = @(x)((L1*x).*(L2*x) - x + F);  % Lorenz '96 dynamical system
gradient_fun = @(x)L96_gradient(x,L1,L2,n);     % Lorenz '96 gradient
x_start = unifrnd(-1,1,n,1);            % random initial condition
L_En4DVar = ACC_Localize(n,r1);         % En4DVar localization matrix for covariance
L_4DVar = ACC_Localize(n,r2);           % 4DVar localization matrix for covariance
L_SqEnKF = ACC_Localize(n,r3);          % SqEnKF localization matrix for covariance
L_EDA = ACC_Localize(n,r4);             % EDA localization matrix for covariance
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
mu = mean(EnsembleSqEnKF,2);
%%

%% experiment + observations + DA
% Run for exp_time, and observe every kth variable at every jump_th 
% model time step. Perform DA with each set of observations.
%%

[EnsembleSqEnKF,EnFEvalSqEnKF] = ODE_RK4_auto_start_Ens(L96fun,EnsembleSqEnKF,dt);
EnsembleEn4DVar = EnsembleSqEnKF(:,1:Ne_En);
X_star_t_En4DVar = mean(EnsembleSqEnKF,2);
X_star_t_EDA = X_star_t_En4DVar;
X_star_t_4DVar = X_star_t_En4DVar;


TimeSeriesSqEnKF = zeros(n,exp_iter);        % array for storing full SqEnKF 
TimeSeries4DVar = zeros(n,exp_iter);         % array for storing full 4DVar
TimeSeriesEDA = zeros(n,exp_iter);           % array for storing full EDA 
TimeSeriesEn4DVar = zeros(n,exp_iter);       % array for storing full En4DVar

spreadVecEn4DVar = spread.*ones(1,exp_iter);   

ErrorVecSqEnKF = zeros(1,exp_iter);
ErrorVec4DVar = zeros(1,exp_iter);
ErrorVecEDA = zeros(1,exp_iter);
ErrorVecEn4DVar = zeros(1,exp_iter);

Time_Series_True = [X,zeros(n,exp_iter)];    % array for storing full true state  
total_steps = 0;

%% from start to first observartions

num_steps = ObsTimes(1)
for ii=2:num_steps
    [Time_Series_True(:,ii),FEvals] = ODE_AB4_auto(Time_Series_True(:,ii-1),FEvals,L96fun,dt);
end

Obs = H*Time_Series_True(:,num_steps) + normrnd(0,sigma,mdim,1);

%% SqEnKF
for jj=2:num_steps
    for mm=1:Ne_Sq
        [EnsembleSqEnKF(:,mm),EnFEvalSqEnKF(:,:,mm)] = ...
            ODE_AB4_auto(EnsembleSqEnKF(:,mm),EnFEvalSqEnKF(:,:,mm),L96fun,dt);
    end
    TimeSeriesSqEnKF(:,jj) = mean(EnsembleSqEnKF,2);
end
[EnsembleSqEnKF,mu_a,~] = DA_SqEnKF(EnsembleSqEnKF,H,Obs,ObsVar,L_SqEnKF,alpha3);
TimeSeriesSqEnKF(:,num_steps) = mu_a;
%%

%% 4DVar
[X_star_t_4DVar,X_star,Time_Series,~,Cov4DVar] = DA_4DVar(X_star_t_4DVar,L96fun,...
    gradient_fun,BCov,H,mu,dt,num_steps-1,Obs,ObsVar,n);
TimeSeries4DVar(:,1:(num_steps-1)) = Time_Series(:,1:(num_steps-1));
%%

%% En4DVar
[EnsembleEn4DVar,X_star_t_En4DVar,~,Time_Series,CovEn4DVar] = DA_En4DVar(X_star_t_En4DVar,...
    EnsembleEn4DVar,L96fun,gradient_fun,H,BCov,mu,dt,jump,Obs,ObsVar,n,L_En4DVar,alpha1);
TimeSeriesEn4DVar(:,1:(num_steps-1)) = Time_Series(:,1:(num_steps-1));
%%

%% EDA
[X_star_t_EDA,spread,CovEDA,Time_Series] = DA_EDA(X_star_t_EDA,L96fun,...
    gradient_fun,BCov,H,mu,dt,jump,Obs,ObsVar,n,L_EDA,alpha4,Ne_EDA);
TimeSeriesEDA(:,1:(num_steps-1)) = Time_Series(:,1:(num_steps-1));
%%

total_steps = total_steps + num_steps;
%% 

%% loop for rest of experiment

for kk=2:q    
    num_steps = ObsTimes(kk)-ObsTimes(kk-1)
    
    for ii=(ObsTimes(kk-1)+1):ObsTimes(kk)
        [Time_Series_True(:,ii),FEvals] = ODE_AB4_auto(Time_Series_True(:,ii-1),FEvals,L96fun,dt);
    end
    
    Obs = H*Time_Series_True(:,ObsTimes(kk)) + normrnd(0,sigma,mdim,1);
    
    %% SqEnKF
    for jj=1:4
        EnsembleSqEnKF = ODE_RK4_auto(EnsembleSqEnKF,L96fun,dt);
        EnFEvalSqEnKF(:,jj,:) = L96fun(EnsembleSqEnKF);
        TimeSeriesSqEnKF(:,ObsTimes(kk-1)+jj) = mean(EnsembleSqEnKF,2);
    end
    
    for jj=5:num_steps
        for mm=1:Ne_Sq
            [EnsembleSqEnKF(:,mm),EnFEvalSqEnKF(:,:,mm)] = ...
                ODE_AB4_auto(EnsembleSqEnKF(:,mm),EnFEvalSqEnKF(:,:,mm),L96fun,dt);
        end
        TimeSeriesSqEnKF(:,ObsTimes(kk-1)+jj) = mean(EnsembleSqEnKF,2);
    end
    [EnsembleSqEnKF,mu_a,~] = DA_SqEnKF(EnsembleSqEnKF,H,Obs,ObsVar,L_SqEnKF,alpha3);
    TimeSeriesSqEnKF(:,ObsTimes(kk)) = mu_a;
    %%
    
    %% 4DVar
    [X_star_t_4DVar,X_star,Time_Series,~,Cov4DVar] = DA_4DVar(X_star_t_4DVar,L96fun,...
        gradient_fun,Cov4DVar,H,X_star_t_4DVar,dt,num_steps,Obs,ObsVar,n);
    TimeSeries4DVar(:,ObsTimes(kk-1):(ObsTimes(kk)-1)) = Time_Series(:,1:(num_steps));
    Cov4DVar = beta*(1+alpha2)*L_4DVar.*Cov4DVar + (1-beta)*BCov;
    Cov4DVar = 0.5*(Cov4DVar + Cov4DVar');
    %%
    
    %% En4DVar
    [EnsembleEn4DVar,X_star_t_En4DVar,spread,Time_Series,CovEn4DVar] = DA_En4DVar(X_star_t_En4DVar,...
        EnsembleEn4DVar,L96fun,gradient_fun,H,CovEn4DVar,mu,dt,jump,Obs,ObsVar,n,L_En4DVar,alpha1);
    TimeSeriesEn4DVar(:,ObsTimes(kk-1):(ObsTimes(kk)-1)) = Time_Series(:,1:(num_steps));
    spreadVecEn4DVar(ObsTimes(kk-1):(ObsTimes(kk)-1)) = spread.*ones(num_steps,1);
    %%
    
    %% EDA
    [X_star_t_EDA,spread,CovEDA,Time_Series] = DA_EDA(X_star_t_EDA,L96fun,...
        gradient_fun,CovEDA,H,X_star_t_EDA,dt,jump,Obs,ObsVar,n,L_EDA,alpha4,Ne_EDA);
    TimeSeriesEDA(:,ObsTimes(kk-1):(ObsTimes(kk)-1)) = Time_Series(:,1:(num_steps));
    %%
end
%%

TimeSeries4DVar(:,end) = X_star_t_4DVar;
TimeSeriesEn4DVar(:,end) = X_star_t_En4DVar;
TimeSeriesEDA(:,end) = X_star_t_EDA;

ErrorSqEnKF = TimeSeriesSqEnKF - Time_Series_True;
Error4DVar = TimeSeries4DVar - Time_Series_True;
ErrorEn4DVar = TimeSeriesEn4DVar - Time_Series_True;
ErrorEDA = TimeSeriesEDA - Time_Series_True;

ErrorVecSqEnKF = vecnorm(ErrorSqEnKF,2)./sqn;
ErrorVec4DVar = vecnorm(Error4DVar,2)./sqn;
ErrorVecEn4DVar = vecnorm(ErrorEn4DVar,2)./sqn;
ErrorVecEDA = vecnorm(ErrorEDA,2)./sqn;


error_parameter_1 = mean(ErrorVec4DVar(10*jump:end));
fprintf('Average RMSE for 4DVar: %g\n',error_parameter_1)
error_parameter_2 = mean(ErrorVecSqEnKF(10*jump:end));
fprintf('Average RMSE for SqEnKF: %g\n',error_parameter_2)
error_parameter_3 = mean(ErrorVecEDA(10*jump:end));
fprintf('Average RMSE for EDA: %g\n',error_parameter_3)
error_parameter_4 = mean(ErrorVecEn4DVar(10*jump:end));
fprintf('Average RMSE for En4DVar: %g\n',error_parameter_4)
fprintf('Average spread: %g\n',mean(spreadVecEn4DVar(10:end)))
%% error plot
set(gcf, 'Position', [25, 25, 1600, 900])
h1 = plot(ErrorVec4DVar,'Color',Color(:,color1),'LineWidth',2.2);
hold on
h2 = plot(ErrorVecSqEnKF,'Color',Color(:,color2),'LineWidth',2.2);
h3 = plot(ErrorVecEDA,'Color',Color(:,color3),'LineWidth',2.2);
h4 = plot(ErrorVecEn4DVar,'Color',Color(:,color4),'LineWidth',2.2);
h5 = plot(spreadVecEn4DVar,'Color',Color(:,color5),'LineWidth',2.2);
title('RMSE & spread')
xlabel('time')
ylabel('RMSE')
legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'4DVar','SqEnKF','EDA','En4DVar','spread')
print('Test_L96_Comparison','-djpeg')
hold off
%%

toc()
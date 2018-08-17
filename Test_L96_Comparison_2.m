tic()           % timer
clc
clear
close all

%% preliminaries
ACC_Colors
warning('off','all')
warning
n = 50;             % dimension of L96 system
sqn = sqrt(n);
Ne_Sq = 40;         % ensemble size
Ne_En = 40;         % ensemble size
Ne_EDA = 40;        % ensemble size
spinup_time = 100;  % for getting onto attractor
exp_time = 40;      % dimensionless time units of DA experiment
long_time = 1000;   % long simulation for creating initial ensemble
dt = 0.02;          % model time step
jump = 10;          % number of model time steps between observations
k = 2;              % observe every kth state variable
F = 8*ones(n,1);    % free parameter on L96 RHS (F = 8 leads to chaotic solutions)
r1 = 7;            % En4DVar localization radius
r2 = 6;             % 4DVar localization radius
r3 = 5.5;           % SqEnKF localization radius
r4 = 7;             % EDA localization radius
alpha1 = 0.15;      % En4DVar inflation parameter
alpha2 = 0.30;      % 4DVar inflation parameter
alpha3 = 0.10;      % SqEnKF inflation parameter
alpha4 = 0.30;      % EDA inflation parameter
ObsVar = 1;         % measurement/observation variance
sigma = sqrt(ObsVar);
beta = 0.1;
color1 = 11;        % SqEnKF
color2 = 29;        % 4DVar
color3 = 3;         % EDA
color4 = 39;        % EnVDvar
color5 = 40;        % True
color6 = 4;         % Obs
color7 = 4;         % spread
color8 = 8;
spinup_iter = floor(spinup_time/dt);    % number of spinup model time steps
exp_iter = floor(exp_time/dt);          % number of experiment model time steps
q = floor(exp_iter/jump);               % number of observed time steps
q_split = floor((3/4)*q);
ObsTimes = jump:jump:(exp_iter+jump);   % vector of times when observation occurs
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
%%

%% experiment + observations + DA
% Run for exp_time, and observe every kth variable at every jump_th 
% model time step. Perform DA with each set of observations.
%%

[EnsembleSqEnKF,EnFEvalSqEnKF] = ODE_RK4_auto_start_Ens(L96fun,EnsembleSqEnKF,dt);

TimeSeriesSqEnKF = zeros(n,exp_iter);        % array for storing full SqEnKF 
TimeSeries4DVar = zeros(n,exp_iter);         % array for storing full 4DVar
TimeSeriesEDA = zeros(n,exp_iter);           % array for storing full EDA 
TimeSeriesEn4DVar = zeros(n,exp_iter);       % array for storing full En4DVar

spreadVecSqEnKF = spread.*ones(1,exp_iter);   

Time_Series_True = [X,zeros(n,exp_iter-1)];     % array for storing full true state  
TimeSeriesObs = zeros(size(H,1),q);             % array for storing all obs
obscounter = 1;
total_steps = 0;

%% from start to first observartions

num_steps = ObsTimes(1);
for ii=2:num_steps
    [Time_Series_True(:,ii),FEvals] = ODE_AB4_auto(Time_Series_True(:,ii-1),FEvals,L96fun,dt);
end

Obs = H*Time_Series_True(:,num_steps) + normrnd(0,sigma,mdim,1);
TimeSeriesObs(:,obscounter) = Obs;
obscounter = obscounter + 1;
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

total_steps = total_steps + num_steps;
%% 

%% loop for spinup of SqEnKF

for kk=2:q_split 
    num_steps = ObsTimes(kk)-ObsTimes(kk-1);
    
    for ii=(ObsTimes(kk-1)+1):ObsTimes(kk)
        [Time_Series_True(:,ii),FEvals] = ODE_AB4_auto(Time_Series_True(:,ii-1),FEvals,L96fun,dt);
    end
    
    Obs = H*Time_Series_True(:,ObsTimes(kk)) + normrnd(0,sigma,mdim,1);
    TimeSeriesObs(:,obscounter) = Obs;
    obscounter = obscounter + 1;
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
    [EnsembleSqEnKF,mu_a,spread,P_a] = DA_SqEnKF(EnsembleSqEnKF,H,Obs,ObsVar,L_SqEnKF,alpha3);
    spreadVecSqEnKF(ObsTimes(kk-1):(ObsTimes(kk)-1)) = spread.*ones(num_steps,1);
    TimeSeriesSqEnKF(:,ObsTimes(kk)) = mu_a;
    %%
end

TimeSeries4DVar(:,1:ObsTimes(q_split)-1) = TimeSeriesSqEnKF(:,1:ObsTimes(q_split)-1);
TimeSeriesEDA(:,1:ObsTimes(q_split)-1) = TimeSeriesSqEnKF(:,1:ObsTimes(q_split)-1);
TimeSeriesEn4DVar(:,1:ObsTimes(q_split)-1) = TimeSeriesSqEnKF(:,1:ObsTimes(q_split)-1);
Cov4DVar = P_a;
CovEDA = P_a;
CovEn4DVar = P_a;
EnsembleEn4DVar = EnsembleSqEnKF(:,1:Ne_En);
X_star_t_En4DVar = mean(EnsembleSqEnKF,2);
X_star_t_EDA = X_star_t_En4DVar;
X_star_t_4DVar = X_star_t_En4DVar;

%% loop for rest of experiment

for kk=q_split+1:q
    num_steps = ObsTimes(kk)-ObsTimes(kk-1);
    
    for ii=(ObsTimes(kk-1)+1):ObsTimes(kk)
        [Time_Series_True(:,ii),FEvals] = ODE_AB4_auto(Time_Series_True(:,ii-1),FEvals,L96fun,dt);
    end
    
    Obs = H*Time_Series_True(:,ObsTimes(kk)) + normrnd(0,sigma,mdim,1);
    TimeSeriesObs(:,obscounter) = Obs;
    obscounter = obscounter + 1;
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
    [EnsembleSqEnKF,mu_a,spread] = DA_SqEnKF(EnsembleSqEnKF,H,Obs,ObsVar,L_SqEnKF,alpha3);
    spreadVecSqEnKF(ObsTimes(kk-1):(ObsTimes(kk)-1)) = spread.*ones(num_steps,1);
    TimeSeriesSqEnKF(:,ObsTimes(kk)) = mu_a;
    %%
    
    %% 4DVar
    [X_star_t_4DVar,X_star,Time_Series,~,Cov4DVar] = DA_4DVar(X_star_t_4DVar,L96fun,...
        gradient_fun,Cov4DVar,H,X_star_t_4DVar,dt,num_steps,Obs,ObsVar,n);
    TimeSeries4DVar(:,ObsTimes(kk-1):(ObsTimes(kk)-1)) = Time_Series(:,1:(num_steps));
    Cov4DVar = (1+alpha2)*L_4DVar.*Cov4DVar;
    Cov4DVar = 0.5*(Cov4DVar + Cov4DVar');
    %%
    
    %% En4DVar
    [EnsembleEn4DVar,X_star_t_En4DVar,~,Time_Series,CovEn4DVar] = DA_En4DVar(X_star_t_En4DVar,...
        EnsembleEn4DVar,L96fun,gradient_fun,H,CovEn4DVar,X_star_t_En4DVar,dt,jump,Obs,ObsVar,n,L_En4DVar,alpha1);
    TimeSeriesEn4DVar(:,ObsTimes(kk-1):(ObsTimes(kk)-1)) = Time_Series(:,1:(num_steps));
    %%
    
    %% EDA
    [X_star_t_EDA,spread,CovEDA,Time_Series] = DA_EDA(X_star_t_EDA,L96fun,...
        gradient_fun,CovEDA,H,X_star_t_EDA,dt,jump,Obs,ObsVar,n,L_EDA,alpha4,Ne_EDA);
    TimeSeriesEDA(:,ObsTimes(kk-1):(ObsTimes(kk)-1)) = Time_Series(:,1:(num_steps));
    %%
end
%%

spreadVecSqEnKF(end) = spread;
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


error_parameter_1 = mean(ErrorVec4DVar(ObsTimes(q_split-1):end));
fprintf('Average RMSE for 4DVar: %g\n',error_parameter_1)
error_parameter_2 = mean(ErrorVecSqEnKF(floor(exp_iter/2):end));
fprintf('Average RMSE for SqEnKF: %g\n',error_parameter_2)
error_parameter_3 = mean(ErrorVecEDA(ObsTimes(q_split-1):end));
fprintf('Average RMSE for EDA: %g\n',error_parameter_3)
error_parameter_4 = mean(ErrorVecEn4DVar(ObsTimes(q_split-1):end));
fprintf('Average RMSE for En4DVar: %g\n',error_parameter_4)
fprintf('Average spread: %g\n',mean(spreadVecSqEnKF(10:end)))
%% error plot 1
figure(1)
set(gcf, 'Position', [25, 25, 1600, 900])
h1 = plot(ErrorVec4DVar,'Color',Color(:,color2),'LineWidth',2.2);
hold on
h2 = plot(ErrorVecEDA,'Color',Color(:,color3),'LineWidth',2.2);
h3 = plot(ErrorVecEn4DVar,'Color',Color(:,color4),'LineWidth',2.2);
h4 = plot(ErrorVecSqEnKF,'Color',Color(:,color1),'LineWidth',2.2);
h5 = plot(spreadVecSqEnKF,'Color',Color(:,color7),'LineWidth',2.2);
title('RMSE & spread')
xlabel('time')
ylabel('RMSE')
legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'4DVar','EDA','En4DVar','SqEnKF','spread')
print('Test_L96_Comparison_1','-djpeg')
hold off
%%

%% error plot 2
xaxis = ObsTimes(q_split-1):exp_iter;
figure(2)
set(gcf, 'Position', [25, 25, 1600, 900])
h1 = plot(xaxis,ErrorVec4DVar(ObsTimes(q_split-1):end),'Color',Color(:,color2),'LineWidth',2.2);
hold on
h2 = plot(xaxis,ErrorVecEDA(ObsTimes(q_split-1):end),'Color',Color(:,color3),'LineWidth',2.2);
h3 = plot(xaxis,ErrorVecEn4DVar(ObsTimes(q_split-1):end),'Color',Color(:,color4),'LineWidth',2.2);
h4 = plot(xaxis,ErrorVecSqEnKF(ObsTimes(q_split-1):end),'Color',Color(:,color1),'LineWidth',2.2);
h5 = plot(xaxis,spreadVecSqEnKF(ObsTimes(q_split-1):end),'Color',Color(:,color7),'LineWidth',2.2);
title('RMSE & spread')
xlabel('time')
ylabel('RMSE')
legend([h1(1),h2(1),h3(1),h4(1),h5(1)],'4DVar','EDA','En4DVar','SqEnKF','spread')
print('Test_L96_Comparison_2','-djpeg')
hold off
%%

%% movie plot 1
Array_SqEnKF = TimeSeriesSqEnKF(:,ObsTimes(q_split-1):end);
Array_4DVar = TimeSeries4DVar(:,ObsTimes(q_split-1):end);
Array_EDA = TimeSeriesEDA(:,ObsTimes(q_split-1):end);
Array_En4DVar = TimeSeriesEn4DVar(:,ObsTimes(q_split-1):end);
Array_True = Time_Series_True(:,ObsTimes(q_split-1):end);
Array_Obs = TimeSeriesObs(:,q_split-1:end);


dim1 = 1;
dim2 = 3;
dim3 = 5;
xx = -7.5;
yy = 11.5;

coords = [xx yy xx yy xx yy];

% L96_movie_1(Array_SqEnKF,Array_4DVar,Array_EDA,Array_En4DVar,Array_True,Array_Obs,...
%     color1,color2,color3,color4,color5,color6,dim1,dim2,dim3,coords,jump)

% L96_movie_2(Time_Series_True,color8,dim1,dim2,dim3,coords)

L96_movie_3(Array_SqEnKF,Array_4DVar,Array_EDA,Array_En4DVar,Array_True,Array_Obs,...
    color1,color2,color3,color4,color5,color6,jump,ObsTimes)
%%

toc()
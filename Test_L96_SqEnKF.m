tic()           % timer
clc
clear
close all

%% preliminaries
ACC_Colors
n = 40;             % dimension of L96 system
Ne = 40;            % ensemble size
spinup_time = 100;  % for getting onto attractor
exp_time = 50;      % dimensionless time units of DA experiment
long_time = 1000;   % long simulation for creating initial ensemble
dt = 0.01;          % model time step
jump = 10;          % number of model time steps between observations
k = 2;              % observe every kth state variable
F = 8*ones(n,1);    % free parameter on L96 RHS (F = 8 leads to chaotic solutions)
r = 4;              % localization radius
alpha = 0.35;       % ensemble inflation parameter
ObsVar = 1;         % measurement/observation variance
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

%% Make ensemble
EnsembleSqEnKF = L96_make_ensemble(L96fun,Ne,dt,long_time,n);
spread = sqrt(trace(cov(EnsembleSqEnKF'))/n);
%%

%% experiment + observations + DA
% Run for exp_time, and observe every kth variable at every jump_th 
% model time step. Perform DA with each set of observations.
%%

[EnsembleSqEnKF,EnFEvalSqEnKF] = ODE_RK4_auto_start_Ens(L96fun,EnsembleSqEnKF,dt);
ErrorVecSqEnKF = zeros(1,exp_iter);
spreadVecSqEnKF = spread*ones(1,exp_iter);
spread = zeros(1,q);
counter = 1;
index = 1;

for kk=1:exp_iter
    [X,FEvals] = ODE_AB4_auto(X,FEvals,L96fun,dt);
    if (counter>1)&&(kk<ObsTimes(counter-1)+5)
        % Right after observations and a DA step, we can't use
        % Adams-Bashforth because our trajectory has been altered.
        % So we restart with RK4 for the first four steps after
        % each DA step.
        EnsembleSqEnKF = ODE_RK4_auto(EnsembleSqEnKF,L96fun,dt);
        EnFEvalSqEnKF(:,index,:) = L96fun(EnsembleSqEnKF);
        index = index + 1;
        fprintf('Performing Runge-Kutta\n')
    else
        % We use Adams-Bashforth to move the model forward until 
        % we get some observations.
        for mm=1:Ne
            [EnsembleSqEnKF(:,mm),EnFEvalSqEnKF(:,:,mm)] = ...
                ODE_AB4_auto(EnsembleSqEnKF(:,mm),EnFEvalSqEnKF(:,:,mm),L96fun,dt);
        end
        index = index + 1;
        fprintf('Performing Adams-Bashforth\n')
    end
    mu_a = mean(EnsembleSqEnKF,2);
    ErrorVecSqEnKF(kk) = norm(mu_a-X,2);
    if kk==ObsTimes(counter)
        Obs = H*X;
        [EnsembleSqEnKF,mu_a,spread] = DA_SqEnKF(EnsembleSqEnKF,H,Obs,ObsVar,L,alpha);
        ErrorVecSqEnKF(kk) = norm(mu_a-X,2);
        spreadVecSqEnKF(ObsTimes(counter):ObsTimes(counter+1)-1) = ...
            spread*ones(1,ObsTimes(counter+1)-ObsTimes(counter));
        counter = counter + 1;
        index = 1;
    end
end

error_parameter = mean(ErrorVecSqEnKF(floor(exp_iter/2):end));
fprintf('Average RMSE: %g\n',error_parameter)
fprintf('Average spread: %g\n',mean(spreadVecSqEnKF(floor(q/2):end)))

%% error plot
set(gcf, 'Position', [25, 25, 1600, 900])
h1 = plot(spreadVecSqEnKF,'Color',Color(:,color1),'LineWidth',2.2);
hold on
h2 = plot(ErrorVecSqEnKF,'Color',Color(:,color2),'LineWidth',2.2);
title('RMSE & spread for SqEnKF')
axis([0 exp_iter 0 3*error_parameter])
xlabel('time')
legend([h1(1),h2(1)],'spread','RMSE')
print('Test_L96_SqEnKF','-djpeg')
hold off
%%

toc()
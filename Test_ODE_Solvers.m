


tic()           % timer
clc
clear
close all

%% preliminaries
ACC_Colors
n = 40;                 % dimension of L96 system
sqn = sqrt(n);
spinup_iter = 10000;    
exp_iter = 100;        % dimensionless time units of DA experiment
dt = 0.01;              % model time step
color1 = 4;
F = 8;
%%

%% setup & utilities
[L1,L2] = L96_get_matrices(n);          % makes matrices for matrix-vector execution of L96
L96fun = @(x)((L1*x).*(L2*x) - x + F);  % Lorenz '96 dynamical system
x_start = unifrnd(-1,1,n,1);            % random initial condition
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

XRK = [X,zeros(n,exp_iter)];
XAB = [X,zeros(n,exp_iter)];
%%

%% experiment + observations + DA
% Run for exp_time, and observe every kth variable at every jump_th 
% model time step. Perform DA with each set of observations.
%%

for kk=1:exp_iter
    [XAB(:,kk+1),FEvals] = ODE_AB4_auto(XAB(:,kk),FEvals,L96fun,dt);
    XRK(:,kk+1) = ODE_RK4_auto(XRK(:,kk),L96fun,dt);
end

error_array = XAB - XRK;
error_vec = zeros(1,exp_iter+1);

for ii=1:exp_iter+1
    error_vec(ii) = norm(error_array(:,ii),2)/sqn;
end

error_parameter = mean(error_vec);
fprintf('Average RMSE: %g\n',error_parameter)

%% error plot
set(gcf, 'Position', [25, 25, 1600, 900])
h1 = plot(error_vec,'Color',Color(:,color1),'LineWidth',2.5);
title('RMSE between AB4 and RK4')
xlabel('time')
ylabel('RMSE')
print('Test_ODE_Solvers','-djpeg')

%%

toc()

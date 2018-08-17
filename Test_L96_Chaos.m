tic()           % timer
clc
clear
close all

%% preliminaries
ACC_Colors
n = 60;                 % dimension of L96 system
sqn = sqrt(n);
spinup_iter = 10000;    
exp_iter = 2000;       % dimensionless time units of DA experiment
dt = 0.02;              % model time step
color1 = 1;
color2 = 13;
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
[X,FEvals1] = ODE_RK4_auto_start(L96fun,x_start,dt);

for ii=1:spinup_iter
    [X,FEvals1] = ODE_AB4_auto(X,FEvals1,L96fun,dt);
end

X_pert = X + [.000001;zeros(n-1,1)];

X1 = [X,zeros(n,exp_iter)];
X2 = [X_pert,zeros(n,exp_iter)];
FEvals2 = FEvals1;
%%

%% experiment + observations + DA
% Run for exp_time, and observe every kth variable at every jump_th 
% model time step. Perform DA with each set of observations.
%%

for kk=1:exp_iter
    [X1(:,kk+1),FEvals1] = ODE_AB4_auto(X1(:,kk),FEvals1,L96fun,dt);
    [X2(:,kk+1),FEvals2] = ODE_AB4_auto(X2(:,kk),FEvals2,L96fun,dt);
end

%% plot
dim1 = 1;
dim2 = 3;
dim3 = 5;
xx = -7.5;
yy = 11.5;

coords = [xx yy xx yy xx yy];

% L96_movie_4(X1,X2,color1,color2,dim1,dim2,dim3,coords)

L96_movie_5(X1,X2,color1,color2)
%%

toc()
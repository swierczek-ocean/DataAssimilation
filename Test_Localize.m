%% preliminaries
ACC_Colors
n = 40;             % dimension of L96 system
r = 4;              % localization radius
alpha = 0.10;       % ensemble inflation parameter
%%

%% setup & utilities
L = ACC_Localize(n,r);                  % localization matrix for covariance
%%

%% test
A = ones(n,n);
surf(L.*A)
%%

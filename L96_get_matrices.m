function [L1,L2] = L96_get_matrices(n)
%% Makes matrices for matrix-vector execution of L96 RHS
% https://en.wikipedia.org/wiki/Lorenz_96_model
% L96: dx(i)/dt = [x(i+1)-x(i-2)]*x(i-1) - x(i) + F
% We wish to implement this as 
% dX/dt = L1*X.*L2*X - X + F with state vector X.
%%

v1 = ones(n-2,1);
v2 = ones(n-1,1);
L1 = diag(v2,1) - diag(v1,-2);
L1(n,1) = 1; L1(1,n-1) = -1; L1(2,n) = -1;
L2 = diag(v2,-1);
L2(1,n) = 1;


end


function gradient = L96_gradient(X,L1,L2,n)
%% Calculation of gradient of L96 RHS
% Don't know what else to say about it.
% The end.
%%

gradient = diag(L1*X)*L2 + diag(L2*X)*L1 - eye(n);

end


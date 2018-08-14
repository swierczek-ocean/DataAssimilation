function [H,m] = L96_get_H(n,k)
%% creates observation operator H
% Works under the restrictive framework of linear 
% observation operator.
% For system of dimension n, observe every kth state
% variable. Works for k = 2, not sure about other values.
%% 

m = floor(n/k);
H = zeros(m,n);

for ii=1:m
   H(ii,k*ii-1)=1; 
end

end


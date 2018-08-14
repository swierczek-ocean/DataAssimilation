function [Ensemble,mu_a,spread] = DA_SQEnKF(Ensemble,H,Obs,ObsVar,L,alpha)
%% Local Ensemble Transform Kalman Filter
% Performs one step of an LETKF data assimilation algorithm.
% 
% The inputs are an ensemble at time k,
% your observations at time k,
% and your model time step dt.
% 
% 
% 'jump' is the number of model time steps between observations
% e.g., if your model has a time step of 0.01 s, and you collect 
% observations every 0.1 s, then jump = 10.
% 
% L and alpha are localization and inflation parameters.
% If you don't know what those are, put alpha = 0 and L = I.
%%
[dim,Ne] = size(Ensemble);
nobs = size(Obs,1);
Rm = sqrt(ObsVar).*eye(nobs);

% for ii=1:jump
%     Ensemble = Model(Ensemble);             % forecast ensemble
% end

P_f = (1+alpha).*L.*cov(Ensemble');         % forecast covariance
mu_f = mean(Ensemble,2);                    % forecast mean
X_f = (Ensemble - mu_f)./sqrt(Ne-1);        % forecast perturbations
K = P_f*H'*((H*P_f*H'+Rm)\eye(nobs));       % Kalman Gain matrix
mu_a = mu_f +K*(Obs-H*mu_f);                % analysis mean
Z = sqrt(1+alpha)*X_f;         
tmp = Z'*H'*(Rm\(H*Z));
tmp = .5*(tmp+tmp');
[E,OM] = eig(tmp);
T = E*(sqrtm(eye(Ne)+OM)\E');               
Za = Z*T;                                   % analysis perturbations
Ensemble = mu_a + sqrt(Ne-1)*Za;            % analysis ensemble
P_a = (eye(n)-K*H)*P_f;                     % analysis covariance
spread = sqrt(trace(P_a)/dim);

end


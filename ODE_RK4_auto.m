function X = ODE_RK4_auto(X,ode_rhs_fun,dt)
%% Fourth order Runge-Kutta method
% Performs one step of RK4
% on ODE: dx/dt = ode_rhs_fun
% with an autonomous RHS
%%

k1=ode_rhs_fun(X);
k2=ode_rhs_fun(X+0.5*dt.*k1);
k3=ode_rhs_fun(X+0.5*dt.*k2);
k4=ode_rhs_fun(X+dt.*k3);
X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);


end


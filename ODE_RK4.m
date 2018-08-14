function X = ODE_RK4(X,ode_rhs_fun,dt,time_start)
%% Fourth order Runge-Kutta method
% Performs one step of RK4
% on ODE: dx/dt = ode_rhs_fun
% with a time dependent RHS
%%

k1=ode_rhs_fun(time_start,X);
k2=ode_rhs_fun(time_start+0.5*dt,X+0.5*dt.*k1);
k3=ode_rhs_fun(time_start+0.5*dt,X+0.5*dt.*k2);
k4=ode_rhs_fun(time_start+dt,X+dt.*k3);
X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);


end


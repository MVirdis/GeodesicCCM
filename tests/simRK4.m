function [ts,xs] = simRK4(dyn, x0, tf, dt, K)
%SIMRK4(dyn, x0, tf, dt, K) Simulate a dynamical system using Runge-Kutta
%4th order discretization method
%   Inputs: dyn => a function that takes as input a state x and
%           (optionally) an input u and returns the derivative xdot
%           x0 => initial state {nx, 1}
%           tf => final time
%           dt => time step used to discretize
%           [K] => optional controller function that takes as input the
%           states and returns a command u. If not specified, u=0
%   Outputs: ts => array of times {1, Nsteps}
%            xs => state trajectory {nx, Nsteps}

if nargin > 4
    nx = size(x0{1},1);
else
    nx = size(x0,1);
end
ts = 0:dt:tf;
N = length(ts);
xs = zeros(nx, N);
if nargin > 4
    xs(:,1) = x0{1};
    xcntr = x0{2};
else
    xs(:,1) = x0;
end

% Runge-Kutta 4th order method
for i = 1:N-1
    x = xs(:,i);

    if nargin > 4
        u_ = K(x,xcntr);
        u = u_{1};
        xcntr = u_{2};
        k1 = dyn(x, u);
        x2 = x + dt*k1/2;
        k2 = dyn(x2, u);
        x3 = x + dt*k2/2;
        k3 = dyn(x3, u);
        x4 = x + dt*k3;
        k4 = dyn(x4, u);
    else
        k1 = dyn(x);
        x2 = x + dt*k1/2;
        k2 = dyn(x2);
        x3 = x + dt*k2/2;
        k3 = dyn(x3);
        x4 = x + dt*k3;
        k4 = dyn(x4);
    end
    
    xs(:,i+1) = xs(:,i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

end

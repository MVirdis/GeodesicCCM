% Filename: Kccm.m
% Author: Mario Virdis
% Created: 2025-04-16
% Description: See function header.
% Version: 1.0

function u = Kccm(config, dK, xeq,x,xk)

u = {0,xk};

% Find geodesic
if config.N > 1
    copt = geodesicFind(config,xeq,x,xk);
else
    copt = geodesicFind(config,xeq,x);
end

% Evaluate along geodesic
u_ = evalDiffK(config,copt,dK);
u = {u_, copt};

end


% Filename: Kccm.m
% Author: Mario Virdis
% Created: 2025-04-16
% Description: See function header.
% Version: 1.0

function u = Kccm(config, dK, xeq,x,xk)

u = {0,xk};

% Find geodesic
params = geodesicParams(config, xeq, x);
opts = optimoptions('fmincon','Algorithm','sqp',...
                    'Display','off', ...
                    'MaxFunctionEvaluations',1e6);
[copt,~,flag] = fmincon(params.f, xk, [], [], ...
                        params.Aeq, params.beq, [], [], [], opts);

% 1  First order optimality conditions satisfied.
% 0  Too many function evaluations or iterations.
% -1  Stopped by output/plot function.
% -2  No feasible point found.
switch(flag)
    case 1
        % Evaluate along geodesic
        u_ = evalDiffK(config,copt,dK);
        u = {u_, copt};
    case 0
        error('Too many function evaluations or iterations.');
    case -1
        error('Stopped by output/plot function.');
    case -2
        error('No feasible point found.');
end

end


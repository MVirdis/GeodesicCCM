function c = geodesicFind(config,xeq,x,c0)
%GEODESICFIND Find a geodesic path

params = geodesicParams(config, xeq, x);

c0_ = params.c0;
if nargin == 4
    c0_ = c0;
end

if config.N > 1 % - Only optimize if the polynomials are at least quadratic
    opts = optimoptions('fmincon','Algorithm','sqp',...
                        'Display','off', ...
                        'MaxFunctionEvaluations',1e6);
    [copt,~,flag] = fmincon(params.f, c0_, [], [], ...
                            params.Aeq, params.beq, [], [], [], opts);
elseif config.N == 1
    % - for N==1 the best geodesic is a straight line connecting x* to x
    copt = c0_;
    flag = 1;
else % - N <= 0 is not a well-posed path in R^n
    error('The degree N must be >= 1')
end

c = c0_;

% 1  First order optimality conditions satisfied.
% 0  Too many function evaluations or iterations.
% -1  Stopped by output/plot function.
% -2  No feasible point found.
switch(flag)
    case 1
        c = copt;
    case 0
        error('Too many function evaluations or iterations.');
    case -1
        error('Stopped by output/plot function.');
    case -2
        error('No feasible point found.');
end

end

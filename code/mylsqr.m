function [x,iter,stats] = mylsqr(matvec,adjvec,b,tol,maxit,varargin)
    summary = [];
    if ~isempty(varargin)
        summary = varargin{1};
    end
    if length(varargin) > 1 && ~isempty(varargin{2}) && norm(varargin{2}) ~= 0
        x = varargin{2};
        b = b - matvec(x);
    else
        x = []; 
    end
    if length(varargin) > 2 && ~isempty(varargin{3})
        verbose = varargin{3};
    else
        verbose = false;
    end
    stats = [];
    bnorm = norm(b);
    beta = bnorm; u = b / beta;
    v = adjvec(u); alpha = norm(v); v = v / alpha;
    w = v;
    if isempty(x)
        x = zeros(size(v));
    end
    phibar = beta; rhobar = alpha;
    resnorm = beta * alpha;
    if ~isempty(summary); stats(end+1,:) = summary(x); end 
    for iter = 1:maxit
        u = matvec(v) - alpha*u; beta = norm(u); u = u / beta;
        v = adjvec(u) - beta*v; alpha = norm(v); v = v / alpha;
        rho = sqrt(rhobar^2 + beta^2);
        c = rhobar / rho;
        s = beta / rho;
        theta = s * alpha;
        rhobar = - c * alpha;
        phi = c * phibar;
        phibar = s * phibar;
        x = x + (phi/rho) * w;
        w = v - (theta/rho) * w;
        if ~isempty(summary); stats(end+1,:) = summary(x); end %#ok<AGROW> 
        if abs(phibar * alpha * c) <= tol * resnorm; break; end
        if verbose
            fprintf('%d\t%e\n',iter,abs(phibar * alpha * c)/resnorm)
        end
    end
end
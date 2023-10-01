function [x,stats] = iterative_sketching(A,b,varargin)
%ITERATIVE_SKETCHING Solve A*x = b in the least-squares sense by the method
%of iterative sketching
%   Optional parameters (use [] for default value):
%   - d: sketching dimension (default max(10*m/n, 20*n) for basic method
%     and max(10*m/n, 4*n) with optimal damping)
%   - q: number of steps. If unspecified, number of steps will be
%     adaptively
%   - summary: a function of the current iterate to be recorded at each
%     iteration. All summary values will be rows of stats, the second
%     output of this function.
%   - verbose: if true, print at each iteration.
%   - damping: damping coefficient (default 0). If 'optimal', the optimal
%     coefficient will be computed as a function of d and n.
%   - momentum: momentum coefficient (default 0). If 'optimal', the optimal
%     coefficient will be computed as a function of d and n.

    m = size(A,1);
    n = size(A,2);
    
    if length(varargin) >= 1 && ~isempty(varargin{1})
        d = varargin{1};
    else
        if length(varargin) >= 5 && ~isempty(varargin{5}) && strcmp(varargin{5}, 'optimal') 
            if length(varargin) >= 6 && ~isempty(varargin{6}) && strcmp(varargin{6}, 'optimal') 
                C = 1;
            else
                C = sqrt(2);
            end
            min_ratio = 4;
        else
            C = 2 + sqrt(2);
            min_ratio = 20;
        end
        d = max(ceil(C^2 * n * exp(lambertw(4*m/n^2*log(1/eps)/C^2))),min_ratio*n);
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        q = varargin{2};
    else
        q = [];
    end
    
    if length(varargin) >= 3 && ~isempty(varargin{3})
        summary = varargin{3};
    else
        summary = [];
    end

    if length(varargin) >= 4 && ~isempty(varargin{4})
        verbose = varargin{4};
    else
        verbose = false;
    end

    if length(varargin) >= 5 && ~isempty(varargin{5})
        damping = varargin{5};
    else
        damping = 1;
    end

    if length(varargin) >= 6 && ~isempty(varargin{6})
        momentum = varargin{6};
    else
        momentum = 0;
    end

    if strcmp(damping, 'optimal') 
        if strcmp(momentum, 'optimal')
            momentum = n/d;
            damping = (1 - momentum)^2;
        elseif momentum == 0
            r = n/d;
            damping = (1-r)^2/(1+r);
        else
            error('Optimal damping for nonzero momentum is not implemented')
        end
    end
    
    stats = [];
    
    S = sparse_sign(d, m, min(8,d));
    B = S*A;
    [Q,R] = qr(B,'econ');
    x = R\(Q'*(S*b));
    xold = x;
    rold = b-A*x;
    if ~isempty(summary); stats(end+1,:) = summary(x); end
    
    Acond = condest(R);
    if isempty(q)
        z = randn(n,1);
        for i = 1:ceil(log(n))
            z = z / norm(z); z = R'*z;
            z = z / norm(z); z = R*z;
        end
        Anorm = norm(z);
    end

    resest = norm(rold) / norm(b);
    if Acond >= 5e-3/eps && resest >= sqrt(eps)
        warning('Condition number (est = %e) and relative residual (est = %e) both appear to be large',...
            Acond, resest)
    end

    iter = 1;
    while true
        xcopy = x;
        x = x + damping * (R\(R'\(A'*rold))) + momentum*(x-xold);
        xold = xcopy;
        r = b-A*x;
        updatenorm = norm(r-rold);
        rold = r;
        if ~isempty(summary); stats(end+1,:) = summary(x); end
        if verbose; fprintf('%d\t%e\n',iter,updatenorm); end
        iter = iter + 1;
        if (~isempty(q) && iter > q) || (isempty(q) &&...
                updatenorm <= eps*(Anorm * norm(x) + 0.01*Acond*norm(r)))
            break
        elseif isempty(q) && iter >= 100
            warning('Iterative sketching failed to meet tolerance')
            break
        end
    end
end
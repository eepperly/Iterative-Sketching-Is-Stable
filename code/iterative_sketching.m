function [x,stats] = iterative_sketching(A,b,varargin)
    m = size(A,1);
    n = size(A,2);
    
    if length(varargin) >= 1 && ~isempty(varargin{1})
        d = varargin{1};
    else
        if ~isempty(varargin{5}) && strcmp(varargin{5}, 'optimal') 
            d = max(ceil(10*m/n), 4*n);
        else
            d = max(ceil(10*m/n), 20*n);
        end
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
    oldupdatenorms = [Inf Inf];
    
    iter = 1;
    while true
        xcopy = x;
        x = x + damping * (R\(R'\(A'*rold))) + momentum*(x-xold);
        xold = xcopy;
        r = b-A*x;
        updatenorm = norm(r - rold);
        rold = r;
        if ~isempty(summary); stats(end+1,:) = summary(x); end
        if verbose; fprintf('%d\t%e\n',iter,updatenorm); end
        iter = iter + 1;
        if (~isempty(q) && iter > q) || (isempty(q) && updatenorm > max(oldupdatenorms) && iter >= 10)
            break
        elseif isempty(q) && iter >= 100
            warning('Iterative sketching failed to meet tolerance')
            break
        end
        oldupdatenorms(2) = oldupdatenorms(1);
        oldupdatenorms(1) = updatenorm;
    end
end
function [x,stats] = bad_iterative_sketching3(A,b,varargin)
    m = size(A,1);
    n = size(A,2);
    
    if length(varargin) >= 1 && ~isempty(varargin{1})
        d = varargin{1};
    else
        d = max(ceil(10*m/n), 20*n);
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
    
    stats = [];
    
    S = sparse_sign(d, m, min(8,d));
    B = S*A;
    [~,R] = qr(B,'econ');
    x = zeros(n,1);
    if ~isempty(summary); stats(end+1,:) = summary(x); end
    rnorm = norm(b-A*x);
    bnorm = norm(b);
    
    iter = 1;
    while true
        x = x + R\(R'\(A'*(b-A*x)));
        rnormnew = norm(b-A*x);
        if ~isempty(summary); stats(end+1,:) = summary(x); end
        fprintf('%d\t%e\n',iter,rnormnew/bnorm);
        iter = iter + 1;
        if (~isempty(q) && iter > q) || (isempty(q) && abs(rnorm-rnormnew) <= 0.5 * eps() * bnorm)
            break
        elseif isempty(q) && iter >= 100
            warning('Iterative sketching failed to meet tolerance')
            break
        end
        rnorm = rnormnew;
    end
end
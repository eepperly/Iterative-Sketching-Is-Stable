function [x,stats] = sketch_and_precondition(A,b,varargin)

    m = size(A,1);
    n = size(A,2);

    if length(varargin) >= 1 && ~isempty(varargin{1})
        d = varargin{1};
    else
        d = 2*n;
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        q = varargin{2};
    else
        q = 100;
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
        opts = varargin{5};
    else
        opts = '';
    end

    S = sparse_sign(d,m,8);
    [Q,R] = qr(S*A,'econ');

    if contains(opts, 'cold')
        y0 = zeros(n,1);
    else
        y0 = Q' * (S*b);
    end

    if ~isempty(summary)
        summary = @(y) summary(R\y);
    end

    if contains(opts, 'cgne')
        [y,~,stats] = mycg(@(y) A'*(A*y),@(y) R\(R'\y),A'*b,0,q,...
            summary,y0,verbose);
    else
        [y,~,stats] = mylsqr(@(y) A*(R\y),@(y) R'\(A'*y),b,0,q,...
            summary,y0,verbose);
    end
    x = R\y;

end
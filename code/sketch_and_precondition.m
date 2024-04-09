function [x,stats] = sketch_and_precondition(A,b,varargin)
%ITERATIVE_SKETCHING Solve A*x = b in the least-squares sense by the
%sketch-and-precondition method
%   Optional parameters (use [] for default value):
%   - d: sketching dimension (default 2*n).
%   - q: number of steps.
%   - summary: a function of the current iterate to be recorded at each
%     iteration. All summary values will be rows of stats, the second
%     output of this function.
%   - verbose: if true, print at each iteration.
%   - opts: specify the solver ('cgne' or 'lsqr') and whether to use a
%     'warm' start (initial iterate given by sketch-and-solve) or a 'cold'
%     start (initial iterate zero). Defaults to 'lsqrwarm'

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

    if length(varargin) >= 6 && ~isempty(varargin{6})
        reproducible = varargin{6};
    else
        reproducible = false;
    end
    
    if ~reproducible && exist('sparsesign','file') == 3
        S = sparsesign(d,m,8);
    else
        warning(['Using slow implementation sparsesign_slow.m. ' ...
            'For the better ' ...
            'implementation, build the mex file `sparsesign.c` ' ...
            'using the command `mex sparsesign.c`.']);
        S = sparsesign_slow(d,m,8);
    end

    [Q,R] = qr(full(S*A),'econ');

    if contains(opts, 'cold')
        y0 = zeros(n,1);
    else
        y0 = Q' * (S*b);
    end

    if contains(opts, 'cgne')
        [x,~,stats] = mycg(@(y) A'*(A*y),@(y) R\(R'\y),A'*b,0,q,...
            summary,R\y0,verbose);
    elseif contains(opts, 'sym')
        if ~isempty(summary)
            summary = @(y) summary(R\y);
        end
        [y,~,stats] = mycg(@(y) R'\(A'*(A*(R\y))),@(y) y,R'\(A'*b),0,q,...
            summary,y0,verbose);
        x = R\y;
    else
        if ~isempty(summary)
            summary = @(y) summary(R\y);
        end
        [y,~,stats] = mylsqr(@(y) A*(R\y),@(y) R'\(A'*y),b,0,q,...
            summary,y0,verbose);
        x = R\y;
    end

end
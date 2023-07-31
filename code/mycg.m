function [x,iter,stats] = mycg(matvec,prec,b,tol,maxit,varargin)
%MYCG Preconditioned conjugate gradient 
% Optional arguments (set to [] for default values):
% 1. summary: function mapping current CG iterate to a row vector of
%    information to be returned in the 'stats' output
% 2. x: initial iterate for CG
% 3. verbose: whether to print convergence statistics.
summary = [];
if ~isempty(varargin)
    summary = varargin{1};
end

if length(varargin) > 1 && ~isempty(varargin{2}) && norm(varargin{2}) ~= 0
    x = varargin{2};
    r = b - matvec(x);
else
    x = zeros(size(b)); 
    r = b;
end

if length(varargin) > 2 && ~isempty(varargin{3})
    verbose = varargin{3};
else
    verbose = false;
end

stats = [];
bnorm = norm(b); rnorm = bnorm;
z = prec(r); p = z;
if ~isempty(summary); stats(end+1,:) = summary(x); end 
for iter = 1:maxit
    if verbose
        fprintf('%d\t%e\n',iter,rnorm/bnorm)
    end
    v = matvec(p);
    zr = z'*r; eta = zr / (v'*p);
    x = x + eta*p;
    r = r - eta*v;
    z = prec(r);
    gamma = z'*r/zr;
    p = z + gamma*p;
    rnorm = norm(r);
    if ~isempty(summary); stats(end+1,:) = summary(x); end %#ok<AGROW> 
    if rnorm <= tol * bnorm; break; end
end
end

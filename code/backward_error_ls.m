function be = backward_error_ls(A, b, y, varargin)
    if isempty(varargin)
        theta = Inf;
    else
        theta = varargin{1};
    end
    
    normy = norm(y);
    
    r = b - A*y;
    normr = norm(r);
    
    if isinf(theta)
        mu = 1;
    else
        mu = theta^2*normy^2; mu = mu / (1 + mu);
    end
    phi = sqrt(mu) * normr / normy;
    
    be = min(phi, min(svd([A phi*(eye(size(A,1))-r*r'/normr^2)])));
end
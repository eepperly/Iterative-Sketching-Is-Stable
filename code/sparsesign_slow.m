function S = sparsesign_slow(d,n,zeta)
    rows = zeros(n*zeta,1);
    for i = 1:n
        rows((i-1)*zeta+1:i*zeta) = randsample(d,zeta);
    end
    cols = kron((1:n)',ones(zeta,1)); % zeta nonzeros per column
    signs = (2*randi(2,n*zeta,1) - 3); % uniform random +/-1 values
    S = sparse(rows, cols, signs / sqrt(zeta), d, n);
end

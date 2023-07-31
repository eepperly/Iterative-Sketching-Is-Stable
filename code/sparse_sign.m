function Phi = sparse_sign(d,N,zeta)
    rows = randi(d,N*zeta,1); % random coordinates in 1,...,d
    cols = kron((1:N)',ones(zeta,1)); % zeta nonzeros per column
    signs = (2*randi(2,N*zeta,1) - 3); % uniform random +/-1 values
    Phi = sparse(rows, cols, signs / sqrt(zeta), d, N);
end
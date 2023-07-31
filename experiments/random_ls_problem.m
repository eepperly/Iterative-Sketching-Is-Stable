function [A,b,x,r] = random_ls_problem(m,n,cond_A,res_size)
U = haarorth(m);
V = haarorth(n);
S = diag(logspace(-log10(cond_A),0,n));
A = U(:,1:n)*S*V';
x = orth(randn(n,1));
r = orth(U(:,(n+1):end)*randn(m-n,1)) * res_size;
b = A*x + r;
end
%% Set up problem

m = 200;
n = 100;

trials = 100; 

cond_A = 10^7;
res_size = 10^-6;

[U,~,V] = svd(randn(m,n));
S = diag(logspace(-log10(cond_A),0,n));
A = U(:,1:n)*S*V';

x = randn(n,1);

r = U(:,(n+1):end)*randn(m-n,1);
r = r / norm(r) * res_size;

b = A*x + r;

%% Sketch and precondition - CG

S = sparse_sign(2*n,m,8);
SA = S*A;
[Q,R] = qr(SA,'econ');
x0 = R\(Q'*(S*b));
[~,~,spre_cg]=mycg(@(y) A'*(A*y),@(y) R\(R'\y),A'*b,0,trials,@(y) [norm(A*y-b)-norm(r);norm(y-x);backward_error_ls(A,b,y)],x0,true);

%% Sketch and precondition - LSQR

[~,~,spre_lsqr]=mylsqr(@(y) A*(R\y), @(y) R'\(A'*y),b,0,trials,@(y) [norm(A*(R\y)-b)-norm(r);norm(R\y-x);backward_error_ls(A,b,R\y)],Q'*(S*b),true);
spre_lsqr = [spre_cg(1,:);spre_lsqr];

%% Sketch and precondition - LSQR from cold start

[~,~,spre_cold]=mylsqr(@(y) A*(R\y), @(y) R'\(A'*y),b,0,trials,@(y) [norm(A*(R\y)-b)-norm(r);norm(R\y-x);backward_error_ls(A,b,R\y)],[],true);
spre_cold = [[norm(b)-norm(r) norm(x) backward_error_ls(A,b,zeros(size(x)))];spre_cold];

%% Iterative sketching

[~,itsk]=iterative_sketching(A,b,20*n,100,@(y) [norm(A*y-b)-norm(r);norm(y-x);backward_error_ls(A,b,y)]);

%% QR

[Q,R] = qr(A,'econ');
y = R\(Q'*b);
qr_vals = [norm(b-A*y)-norm(r);norm(x-y);backward_error_ls(A,b,y)];

%% Plot

close all
for j = 1:3
    figure(j)
    semilogy(0:trials, spre_cg(:,j)); hold on
    semilogy(0:trials, spre_cold(:,j));
    semilogy(0:trials, spre_lsqr(:,j));
    semilogy(0:trials, itsk(:,j));
    yline(qr_vals(j))
    xlabel('Iteration')
    if j == 1
        legend('CGNE (warm)', 'LSQR (cold)', 'LSQR (warm)', 'IS (warm)')
        ylabel('Residual $\|Ax-b\|$')
%         saveas(gcf,'residual.fig')
%         saveas(gcf,'residual.png')
    elseif j == 2
        ylabel('Forward error $\|x-x_\star\|$')
%         saveas(gcf,'forward.fig')
%         saveas(gcf,'forward.png')
    elseif j == 3
        ylabel('Backward error $\eta(x)$')
%         saveas(gcf,'backward.fig')
%         saveas(gcf,'backward.png')
    end
end
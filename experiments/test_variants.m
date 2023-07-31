%% Set up problem

addpath('../code')

m = 4000;
n = 50;

trials = 100; 

cond_A = 10^10;
res_size = 10^-2;

U = haarorth(m);
V = haarorth(n);
S = diag(logspace(-log10(cond_A),0,n));
A = U(:,1:n)*S*V';
x = orth(randn(n,1));

r = orth(U(:,(n+1):end)*randn(m-n,1)) * res_size;

b = A*x + r;

% summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);backward_error_ls(A,b,y)/norm(A,'fro')];
summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);0];

%% Damped

[~,damping]=iterative_sketching(A,b,12*n,trials,summary,true,'optimal');

%% Momentum

[~,momentum]=iterative_sketching(A,b,12*n,trials,summary,true,'optimal','optimal');

%% Iterative sketching

[~,itsk]=iterative_sketching(A,b,12*n,trials,summary,true);

%% QR

y = A\b;
qr_vals = summary(y);

%% Plot

close all
for j = 1:3
    figure(j)
    semilogy(0:trials, itsk(:,j), 'LineWidth', 4);hold on
    semilogy(0:trials, damping(:,j), '-.', 'LineWidth', 4); 
    semilogy(0:trials, momentum(:,j), ':', 'LineWidth', 4); 
    yline(qr_vals(j),'k:', 'LineWidth', 3)
    xlabel('Iteration $i$')
    if j == 1
        ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
        legend('Plain',...
            'Damped',...
            'Momentum',...
            'QR')
        saveas(gcf,'../figs/variants_forward.fig')
        saveas(gcf,'../figs/variants_forward.png')
    elseif j == 2
        ylabel('Residual error $\|\mbox{\boldmath $r$}(\mbox{\boldmath $x$})-\mbox{\boldmath $r$}(\mbox{\boldmath $\widehat{x}$}_i)\|/\|\mbox{\boldmath $b$}\|$')
        saveas(gcf,'../figs/variants_residual.fig')
        saveas(gcf,'../figs/variants_residual.png')
    elseif j == 3
        ylabel('Backward error $\eta(\mbox{\boldmath $\widehat{x}$}_i)//\|\mbox{\boldmath $A$}\|_{\rm F}$')
        saveas(gcf,'../figs/variants_backward.fig')
        saveas(gcf,'../figs/variants_backward.png')
    end
end
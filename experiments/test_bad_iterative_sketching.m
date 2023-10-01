%% Set up problem

addpath("bad_implementations/")
addpath("../code")

m = 4000;
n = 50;

trials = 80; 

cond_A = 10^10;
res_size = 10^-6;
[A,b,x,r] = random_ls_problem(m, n, cond_A, res_size);

real_run = true;
if real_run
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);...
        backward_error_ls(A,b,y)/norm(A,'fro')];
else
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);0];
end

%% Bad approach 1

fprintf('Bad approach 1\n')
[~,bad1]=bad_iterative_sketching1(A,b,20*n,trials,summary);

%% Bad approach 2

fprintf('Bad approach 2\n')
[~,bad2]=bad_iterative_sketching2(A,b,20*n,trials,summary);

%% Bad approach 3

fprintf('Bad approach 3\n')
[~,bad3]=bad_iterative_sketching3(A,b,20*n,trials,summary);

%% Iterative sketching

fprintf('Good approach\n')
[~,itsk]=iterative_sketching(A,b,20*n,trials,summary,true);

%% QR

fprintf('QR\n')
[Q,R] = qr(A, 'econ');
y = R\(Q'*b);
qr_vals = summary(y);

%% Plot

close all
for j = 1:3
    figure(j)
    semilogy(0:trials, itsk(:,j), 'LineWidth', 4);hold on
    semilogy(0:trials, bad2(:,j), '-.', 'LineWidth', 4); 
    semilogy(0:trials, bad1(:,j), ':', 'LineWidth', 4);
    semilogy(0:trials, bad3(:,j), '--', 'LineWidth', 4);
    yline(qr_vals(j),'k:', 'LineWidth', 3)
    xlabel('Iteration $i$')
    if j == 1
        ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
        axis([-Inf Inf 1e-6 1e10])
        legend('Stable (Algorithm 2)',...
            'Bad matrix: $(\mbox{\boldmath $SA$})^\top(\mbox{\boldmath $SA$})$',...
            'Bad residual: $\mbox{\boldmath $A^\top b$}-\mbox{\boldmath $A^\top$} (\mbox{\boldmath $A\widehat{x}$}_i)$',...
            'Bad initialization: $\mbox{\boldmath $x$}_0=\mbox{\boldmath $0$}$','QR')
        if real_run
            saveas(gcf,'../figs/bad_forward.fig')
            saveas(gcf,'../figs/bad_forward.png')
        end
    elseif j == 2
        ylabel('Residual error $\|\mbox{\boldmath $r$}(\mbox{\boldmath $x$})-\mbox{\boldmath $r$}(\mbox{\boldmath $\widehat{x}$}_i)\|/\|\mbox{\boldmath $b$}\|$')
        axis([-Inf Inf 1e-16 1e10])
        if real_run
            saveas(gcf,'../figs/bad_residual.fig')
            saveas(gcf,'../figs/bad_residual.png')
        end
    elseif j == 3
        ylabel('Backward error $\eta(\mbox{\boldmath $\widehat{x}$}_i)//\|\mbox{\boldmath $A$}\|_{\rm F}$')
        if real_run
            saveas(gcf,'../figs/bad_backward.fig')
            saveas(gcf,'../figs/bad_backward.png')
        end
    end
end


%% Save

if real_run
    save('../data/results_bad_iterative_sketching.mat', 'bad1', 'bad2', 'itsk', 'bad3','qr_vals', 'trials')
end


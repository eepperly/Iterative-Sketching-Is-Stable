%% Set up problem

addpath('../code')

m = 4000;
n = 50;

trials = 40; 

cond_A = 10^10;
res_size = 10^-6;

[A,b,x,r] = random_ls_problem(m,n,cond_A,res_size);

real_run = true;
if real_run
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);backward_error_ls(A,b,y)/norm(A,'fro')];
else
    summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);0];
end

%% LSQR Warm

[~,lsqrwarm]=sketch_and_precondition(A,b,20*n,trials,summary,true,'lsqrwarm');

%% LSQR Cold

[~,lsqrcold]=sketch_and_precondition(A,b,20*n,trials,summary,true,'lsqrcold');

%% Iterative Sketching

[~,itsk]=iterative_sketching(A,b,20*n,trials,summary,true);

%% Iterative Sketching with Damping

[~,damp]=iterative_sketching(A,b,20*n,trials,summary,true,'optimal');

%% Iterative Sketching with Momentum

[~,mom]=iterative_sketching(A,b,20*n,trials,summary,true,'optimal','optimal');

%% QR

[Q,R] = qr(A,'econ');
y = R\(Q'*b);
qr_vals = summary(y);

%% Plot

close all
for j = 1:3
    figure(j)
    semilogy(0:trials, lsqrcold(:,j), '-.', 'LineWidth', 4, 'Color', "#D95319"); hold on
    semilogy(0:trials, itsk(:,j), ':', 'LineWidth', 4, 'Color', "#EDB120"); 
    semilogy(0:trials, damp(:,j), '^', 'LineWidth', 1, 'Color',...
        "#7E2F8E", 'MarkerSize', 10, 'MarkerFaceColor', "#7E2F8E"); 
    semilogy(0:trials, mom(:,j), 'o', 'LineWidth', 1, 'Color',...
        "#77AC30",'MarkerSize', 10, 'MarkerFaceColor', "#77AC30"); 
    semilogy(0:trials, lsqrwarm(:,j), 'LineWidth', 4, 'Color', "#0072BD");
    yline(qr_vals(j),'k:', 'LineWidth', 3)
    xlabel('Iteration $i$')
    if j == 1
        ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
        legend('S\&P ($\mbox{\boldmath $x$}_0=\mbox{\boldmath $0$}$)',...
            'IS',...
            'IS+damping',...
            'IS+momentum',...
            'S\&P',...
            'QR')
	if real_run
            saveas(gcf,'../figs/sketch_precondition_forward.fig')
            saveas(gcf,'../figs/sketch_precondition_forward.png')
	end
    elseif j == 2
        ylabel('Residual error $\|\mbox{\boldmath $r$}(\mbox{\boldmath $x$})-\mbox{\boldmath $r$}(\mbox{\boldmath $\widehat{x}$}_i)\|/\|\mbox{\boldmath $b$}\|$')
	if real_run
            saveas(gcf,'../figs/sketch_precondition_residual.fig')
            saveas(gcf,'../figs/sketch_precondition_residual.png')
	end
    elseif j == 3
        ylabel('Backward error $\eta(\mbox{\boldmath $\widehat{x}$}_i)//\|\mbox{\boldmath $A$}\|_{\rm F}$')
	if real_run
            saveas(gcf,'../figs/sketch_precondition_backward.fig')
            saveas(gcf,'../figs/sketch_precondition_backward.png')
	end
    end
end

%% Save

if real_run
    save('../data/results_variants.mat', 'lsqrwarm', 'lsqrcold', 'itsk', 'damp', 'mom', 'qr_vals', 'trials')
end

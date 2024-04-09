%% Set up problem

addpath('../code')

rng(121842)

m = 4000;
n = 50;

trials = 80; 

conds = [1e1 1e10 1e15];
res_sizes = [1e-3 1e-12];

itsk = cell(length(conds), length(res_sizes));
qrs = cell(length(conds), length(res_sizes));
colors = ["#0072BD","#D95319","#EDB120"];
markers = '^so';

real_run = true;

for idx1 = 1:length(conds)
    cond_A = conds(idx1);
    for idx2 = 1:length(res_sizes)
        res_size = res_sizes(idx2);
	    [A,b,x,r] = random_ls_problem(m,n,cond_A,res_size);

        if real_run
            summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);backward_error_ls(A,b,y)/norm(A,'fro')];
        else
            summary = @(y) [norm(y-x)/norm(x);norm(b-A*y-r)/norm(b);0];
        end

        [~,itsk{idx1,idx2}]=iterative_sketching(A,b,20*n,trials,summary,true,[],[],true);
        [Q,R] = qr(A,'econ');
        y = R\(Q'*b);
        qrs{idx1,idx2} = summary(y);
    end
end

%% Plot

close all
for idx2 = 1:length(res_sizes)
    res_size = res_sizes(idx2);
    for idx1 = length(conds):-1:1
        cond_A = conds(idx1);

        for j = 1:3
            figure(3*(idx2-1)+j)
            itsk_idx = itsk{idx1,idx2};
            qr_idx = qrs{idx1,idx2};
            semilogy(0:trials, itsk_idx(1:(trials+1),j), 'LineWidth', 4,...
                'Color', colors{idx1}); hold on
            semilogy(0:trials, itsk_idx(1:(trials+1),j), 'LineWidth', 1,... ...
                'LineStyle', "none",...
                'Color', colors{idx1},'MarkerFaceColor',colors{idx1},...
                'Marker',markers(idx1),'MarkerIndices',1:10:(trials+1),...
                'MarkerSize',15); hold on
            yline(qr_idx(j),':', 'LineWidth', 3, 'Color', colors{idx1})
        end
    end
end
if real_run
    save("../data/results_test_iterative_sketching.mat")
end

%% Save

if real_run
    for idx2 = 1:length(res_sizes)
        figure(3*(idx2-1)+1)
        xlabel('Iteration $i$')
        ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
    
        if idx2 == 1
            legend({'','$\kappa = 10^{15}$','','','$\kappa = 10^{10}$','','','$\kappa = 10^{1}$',''},'Location','best')
        end

	if real_run
	   saveas(gcf, sprintf('../figs/r%d_forward.fig', round(-log10(res_sizes(idx2)))))
           saveas(gcf, sprintf('../figs/r%d_forward.png', round(-log10(res_sizes(idx2)))))
	end
    
        figure(3*(idx2-1)+2)
        xlabel('Iteration $i$')
        ylabel('Residual error $\|\mbox{\boldmath $r$}(\mbox{\boldmath $x$})-\mbox{\boldmath $r$}(\mbox{\boldmath $\widehat{x}$}_i)\|/\|\mbox{\boldmath $b$}\|$')
	if real_run
           saveas(gcf, sprintf('../figs/r%d_residual.fig', round(-log10(res_sizes(idx2)))))
           saveas(gcf, sprintf('../figs/r%d_residual.png', round(-log10(res_sizes(idx2)))))
	end
    
        figure(3*(idx2-1)+3)
        xlabel('Iteration $i$')
        ylabel('Backward error $\mbox{BE}(\mbox{\boldmath $\widehat{x}$}_i)$')
	if real_run
           saveas(gcf, sprintf('../figs/r%d_backward.fig', round(-log10(res_sizes(idx2)))))
           saveas(gcf, sprintf('../figs/r%d_backward.png', round(-log10(res_sizes(idx2)))))
	end
    end
end

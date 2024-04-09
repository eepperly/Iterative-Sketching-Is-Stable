%% Load SUSY data

addpath('../code')
load('../data/susy.mat')
X = X(1:4e6,:); b = b(1:4e6);
bandwidth = 4;
S = randsample(4e6,1e3);
A = exp(-pdist2(X,X(S,:),"euclidean").^2 / (2*bandwidth^2));

sizes = round(logspace(1,3,11));

real_run = true;

%% Run trials

qrs = zeros(length(sizes), 1);
wedins = zeros(length(sizes), 1);
itsks = zeros(length(sizes), 2);
damps = zeros(length(sizes), 2);
moms = zeros(length(sizes), 2);

for i = 1:11
    sizes(i)
    AA = A(:,1:sizes(i));

    % QR
    tic; x_qr = AA\b; qrs(i) = toc;
    qrs(i)

    % Iterative sketching
    tic; x = iterative_sketching(AA,b); itsks(i,1) = toc;
    itsks(i,2) = norm(x - x_qr) / norm(x_qr);
    itsks(i,1)

    % Iterative sketching with damping
    tic; x = iterative_sketching(AA,b,[],[],[],[],'optimal');
    damps(i,1) = toc;
    damps(i,2) = norm(x - x_qr) / norm(x_qr);
    damps(i,1)

    % Iterative sketching with momentum
    tic; x = iterative_sketching(AA,b,[],[],[],[],'optimal','optimal');
    moms(i,1) = toc;
    moms(i,2) = norm(x - x_qr) / norm(x_qr);
    moms(i,1)
    
    % Wedin
    [~,R] = qr(AA,'econ');
    r_qr = b - AA*x_qr;
    cond_A = cond(R)
    norm_A = norm(R)
    wedins(i) = 2.23 * cond_A * (1 + cond_A / norm_A / norm(x_qr) * norm(r_qr)) * eps;
end

if real_run
    save('../data/results_timing.mat', 'qrs', 'itsks', 'damps', 'moms', 'sizes', 'wedins')
end

%% Plot

close all
figure(1)
loglog(sizes, qrs(:,1), 'k', 'LineWidth', 4); hold on
loglog(sizes, itsks(:,1), ':', 'LineWidth',4,'Color',"#EDB120");
xlabel('Number of columns $n$')
ylabel('Time (sec)')
legend({'\texttt{mldivide}', 'Iterative sketching'},'Location','best')

if real_run
    saveas(gcf, '../figs/susy_times_slim.fig')
    saveas(gcf, '../figs/susy_times_slim.png')
end

loglog(sizes, damps(:,1), '^', 'LineWidth',1,'Color',"#7E2F8E",...
    'MarkerFaceColor',"#7E2F8E",'MarkerSize',10);
loglog(sizes, moms(:,1), 'o', 'LineWidth',1,'Color',"#77AC30",...
    'MarkerFaceColor',"#77AC30",'MarkerSize',10);
legend({'QR', 'IS', 'IS+damping', 'IS+momentum'},'Location','best')

if real_run
    saveas(gcf, '../figs/susy_times.fig')
    saveas(gcf, '../figs/susy_times.png')
end

figure(2)
loglog(sizes, wedins, ':', 'LineWidth', 3, 'Color', "#A2142F"); hold on
loglog(sizes, itsks(:,2), ':', 'LineWidth',4,'Color',"#EDB120");
xlabel('Number of columns $n$')
ylabel('Forward error $\|\mbox{\boldmath $\widehat{x}$}_{\mbox{dir}}-\mbox{\boldmath $\widehat{x}$}\|/\|\mbox{\boldmath $\widehat{x}$}_{\mbox{dir}}\|$')

if real_run
    saveas(gcf, '../figs/susy_errors_slim.fig')
    saveas(gcf, '../figs/susy_errors_slim.png')
end

loglog(sizes, damps(:,2), '^', 'LineWidth',1,'Color',"#7E2F8E",...
    'MarkerFaceColor',"#7E2F8E",'MarkerSize',10);
loglog(sizes, moms(:,2), 'o', 'LineWidth',1,'Color',"#77AC30",...
    'MarkerFaceColor',"#77AC30",'MarkerSize',10);

if real_run
    saveas(gcf, '../figs/susy_errors.fig')
    saveas(gcf, '../figs/susy_errors.png')
end
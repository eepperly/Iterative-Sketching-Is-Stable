%% Load sparse data

addpath('../code')
N = 3e6;
b = randn(N,1);
sizes = round(logspace(1,3,11));

real_run = true;

%% Run trials

qrs = zeros(length(sizes), 1);
wedins = zeros(length(sizes), 1);
itsks = zeros(length(sizes), 1);
damps = zeros(length(sizes), 1);
moms = zeros(length(sizes), 1);

for i = 1:length(sizes)
    sizes(i)
    AA = sparsesign(sizes(i),N,3)'*sqrt(3);

    % QR
    fprintf('Direct\n')
    tic; x_qr = AA\b; qrs(i) = toc;
    qrs(i)

    % Iterative sketching
    fprintf('Iterative sketching\n')
    tic; x = iterative_sketching(AA,b); itsks(i,1) = toc;
    itsks(i,1)

    % Iterative sketching with damping
    fprintf('Iterative sketching + damping\n')
    tic; x = iterative_sketching(AA,b,[],[],[],[],'optimal');
    damps(i,1) = toc;
    damps(i,1)

    % Iterative sketching with momentum
    fprintf('Iterative sketching + momentum\n')
    tic; x = iterative_sketching(AA,b,[],[],[],[],'optimal','optimal');
    moms(i,1) = toc;
    moms(i,1)
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
    saveas(gcf, '../figs/sparse_times_slim.fig')
    saveas(gcf, '../figs/sparse_times_slim.png')
end

loglog(sizes, damps(:,1), '^', 'LineWidth',1,'Color',"#7E2F8E",...
    'MarkerFaceColor',"#7E2F8E",'MarkerSize',10);
loglog(sizes, moms(:,1), 'o', 'LineWidth',1,'Color',"#77AC30",...
    'MarkerFaceColor',"#77AC30",'MarkerSize',10);
legend({'\texttt{mldivide}', 'IS', 'IS+damping', 'IS+momentum'},'Location','best')

if real_run
    saveas(gcf, '../figs/sparse_times.fig')
    saveas(gcf, '../figs/sparse_times.png')
end
%% Load SUSY data

load('../data/susy.mat')
bandwidth = 4;
A = exp(-pdist2(X,X(S,:),"euclidean").^2 / (2*bandwidth^2));

sizes = round(logspace(1,3,11));
solutions = cell(11,1);
residuals = cell(11,1);

%% Run trials

qrs = zeros(length(sizes), 3);
itsks = zeros(length(sizes), 3);
damps = zeros(length(sizes), 3);
moms = zeros(length(sizes), 3);

for i = 1:length(sizes)
    sizes(i)
    AA = A(:,1:sizes(i));

    % QR
    tic; x = AA\b; qrs(i,1) = toc;
    qrs(i,2) = norm(x - solutions{i}) / norm(solutions{i});
    qrs(i,3) = norm((b-AA*x)-residuals{i}) / norm(residuals{i});

    % Iterative sketching
    tic; x = iterative_sketching(AA,b); itsks(i,1) = toc;
    itsks(i,2) = norm(x - solutions{i}) / norm(solutions{i});
    itsks(i,3) = norm((b-AA*x)-residuals{i}) / norm(residuals{i});

    % Iterative sketching with damping
    tic; x = iterative_sketching(AA,b,[],[],[],[],'optimal');
    damps(i,1) = toc;
    damps(i,2) = norm(x - solutions{i}) / norm(solutions{i});
    damps(i,3) = norm((b-AA*x)-residuals{i}) / norm(residuals{i});

    % Iterative sketching with momentum
    tic; x = iterative_sketching(AA,b,[],[],[],[],'optimal','optimal');
    moms(i,1) = toc;
    moms(i,2) = norm(x - solutions{i}) / norm(solutions{i});
    moms(i,3) = norm((b-AA*x)-residuals{i}) / norm(residuals{i});
end

%% Plot

close all
figure(1)
loglog(sizes, qrs(:,1), 'LineWidth', 4); hold on
loglog(sizes, itsks(:,1), ':', 'LineWidth',4);
loglog(sizes, damps(:,1), '-.', 'LineWidth',4);
loglog(sizes, moms(:,1), '--', 'LineWidth',4);
% axis([10 1000 0.7 2e3])
xlabel('Number of columns $n$')
ylabel('Time (sec)')
legend('QR', 'IS', 'IS+damping', 'IS+momentum')

figure(2)
loglog(sizes, qrs(:,2), 'LineWidth', 4); hold on
loglog(sizes, itsks(:,2), ':', 'LineWidth',4);
loglog(sizes, damps(:,2), '-.', 'LineWidth',4);
loglog(sizes, moms(:,2), '--', 'LineWidth',4);
xlabel('Number of columns $n$')
ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
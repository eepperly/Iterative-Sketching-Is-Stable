%% Load SUSY data

load('../data/susy.mat')
bandwidth = 4;
A = exp(-pdist2(X,X(S,:),"euclidean").^2 / (2*bandwidth^2));

sizes = round(logspace(1,3,11));

load("../data/true_solutions.mat")
solutions = { x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11 };
residuals = cell(11,1);
for i = 1:11
    residuals{i} = b - A(:,1:sizes(i)) * solutions{i};
end

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
    qrs(i,1)

    % Iterative sketching
    tic; x = iterative_sketching(AA,b); itsks(i,1) = toc;
    itsks(i,2) = norm(x - solutions{i}) / norm(solutions{i});
    itsks(i,3) = norm((b-AA*x)-residuals{i}) / norm(residuals{i});
    itsks(i,1)

    % Iterative sketching with damping
    tic; x = iterative_sketching(AA,b,[],[],[],[],'optimal');
    damps(i,1) = toc;
    damps(i,2) = norm(x - solutions{i}) / norm(solutions{i});
    damps(i,3) = norm((b-AA*x)-residuals{i}) / norm(residuals{i});
    damps(i,1)

    % Iterative sketching with momentum
    tic; x = iterative_sketching(AA,b,[],[],[],[],'optimal','optimal');
    moms(i,1) = toc;
    moms(i,2) = norm(x - solutions{i}) / norm(solutions{i});
    moms(i,3) = norm((b-AA*x)-residuals{i}) / norm(residuals{i});
    moms(i,1)
end

%% Plot

close all
figure(1)
loglog(sizes, qrs(:,1), 'k', 'LineWidth', 4); hold on
loglog(sizes, itsks(:,1), ':', 'LineWidth',4,'Color',"#EDB120");
loglog(sizes, damps(:,1), '^', 'LineWidth',1,'Color',"#7E2F8E",...
    'MarkerFaceColor',"#7E2F8E",'MarkerSize',10);
loglog(sizes, moms(:,1), 'o', 'LineWidth',1,'Color',"#77AC30",...
    'MarkerFaceColor',"#77AC30",'MarkerSize',10);
% axis([10 1000 0.7 2e3])
xlabel('Number of columns $n$')
ylabel('Time (sec)')
legend('QR', 'IS', 'IS+damping', 'IS+momentum')

figure(2)
loglog(sizes, qrs(:,2), 'k', 'LineWidth', 4); hold on
loglog(sizes, itsks(:,2), ':', 'LineWidth',4,'Color',"#EDB120");
loglog(sizes, damps(:,2), '^', 'LineWidth',1,'Color',"#7E2F8E",...
    'MarkerFaceColor',"#7E2F8E",'MarkerSize',10);
loglog(sizes, moms(:,2), 'o', 'LineWidth',1,'Color',"#77AC30",...
    'MarkerFaceColor',"#77AC30",'MarkerSize',10);
xlabel('Number of columns $n$')
ylabel('Forward error $\|\mbox{\boldmath $x$}-\mbox{\boldmath $\widehat{x}$}_i\|/\|\mbox{\boldmath $x$}\|$')
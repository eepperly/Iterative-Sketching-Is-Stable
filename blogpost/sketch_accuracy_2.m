n = 1e5;
k = 5e1;
A = sprand(n,k,1e-2);
Q1 = orth(full(A));
Q2 = orth(randn(n,k));
U1 = orth(randn(50));
U2 = orth(randn(50));
U3 = orth(randn(50));
Q3 = zeros(50^3,k);
for i = 1:k
    Q3(:,i) = kron(U1(:,i),kron(U2(:,i),U3(:,i)));
end
Q4 = [eye(k);zeros(n-k,k)];
Qs = {Q1,Q2,Q3,Q4};
good_parameter_setting = true;

ds = round(logspace(2,4,7));
dis = zeros(length(ds),4);
trials = 100;
addpath('../code')

for d_idx = 1:length(ds)
    d = ds(d_idx)
    for i = 1:4
        Q = Qs{i};
        for trial = 1:trials
            if good_parameter_setting
                S = sparsesign(d,size(Q,1),max(8,ceil(sqrt(4*d/k))));
            else
                S = sparsesign(d,size(Q,1),8);
            end
            SQ = full(S*Q);
            svals = svd(SQ);
            dis(d_idx,i) = dis(d_idx,i)+max(max(svals)-1,1-min(svals))/trials;
        end
    end
end

%% Plot
close all
loglog(ds,dis(:,1),'^:','LineWidth',2,'MarkerFaceColor',"#0072BD",'MarkerSize',16); hold on
loglog(ds,dis(:,2),'s:','LineWidth',2,'MarkerFaceColor',"#D95319",'MarkerSize',16)
loglog(ds,dis(:,3),'o:','LineWidth',2,'MarkerFaceColor',"#EDB120",'MarkerSize',16)
loglog(ds,dis(:,4),'v:','LineWidth',2,'MarkerFaceColor',"#7E2F8E",'MarkerSize',16)
loglog(ds,sqrt(k./ds),'k--','LineWidth',3)
xlabel('Sketching dimension $d$')
ylabel('Distortion $\varepsilon$')
legend({'Sparse', 'Dense', 'Khatri--Rao', 'Identity', '$\sqrt{k/d}$'})
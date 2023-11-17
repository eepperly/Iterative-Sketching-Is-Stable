n = 1e5;
k = 2e2;
Q = [eye(k);zeros(n-k,k)];
trials = 100;

zetas = [1 2 4 8];
ds = 2*round(logspace(2,4,7));
dis = zeros(length(ds),length(zetas));


for d_idx = 1:length(ds)
    d = ds(d_idx)

    for trial = 1:trials
        for zeta_idx = 1:length(zetas)
            zeta = zetas(zeta_idx);
            S = sparsesign(d,n,zeta);
            SQ = S*Q;
            svals = svd(SQ);
            dis(d_idx,zeta_idx) = dis(d_idx,zeta_idx) + max(max(svals)-1,1-min(svals))/trials;
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
legend({'CountSketch', '$\zeta=2$', '$\zeta=4$', '$\zeta=8$', '$\sqrt{k/d}$'},'Location','best')
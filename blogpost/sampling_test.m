n = 1e5;
k = 5e1;
Q = [eye(k);zeros(n-k,k)];
trials = 100;

uniform = zeros(length(ds),1);
leverage = zeros(length(ds),1);
sse = zeros(length(ds),1);

ds = round(logspace(2,4,7));
for d_idx = 1:length(ds)
    d = ds(d_idx)

    for trial = 1:trials
        S = sparsesign(d,n,max(8,ceil(2*sqrt(d/k))));
        SQ = S*Q;
        svals = svd(SQ);
        sse(d_idx) = sse(d_idx) + max(max(svals)-1,1-min(svals))/trials;

        SQ = sqrt(n/d)*Q(randsample(n,d,true),:);
        svals = svd(SQ);
        uniform(d_idx) = uniform(d_idx) + min(max(max(svals)-1,1-min(svals)),1)/trials;

        SQ = sqrt(k/d)*Q(randsample(k,d,true),:);
        svals = svd(SQ);
        leverage(d_idx) = leverage(d_idx) + max(max(svals)-1,1-min(svals))/trials;
    end
end

%% Plot
close all
loglog(ds,uniform,'^:','LineWidth',2,'MarkerFaceColor',"#0072BD",'MarkerSize',16); hold on
loglog(ds,leverage,'s:','LineWidth',2,'MarkerFaceColor',"#D95319",'MarkerSize',16)
loglog(ds,sse,'o:','LineWidth',2,'MarkerFaceColor',"#EDB120",'MarkerSize',16)
loglog(ds,sqrt(k./ds),'k--','LineWidth',3)
xlabel('Sketching dimension $d$')
ylabel('Distortion $\varepsilon$')
legend({'Uniform', 'Leverage scores', 'Sparse sign (suggested $\zeta$)', '$\sqrt{k/d}$'})
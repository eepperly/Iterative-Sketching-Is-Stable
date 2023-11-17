n = 1e5;
k = 5e1;
A = sprand(n,k,1e-2);
Q = orth(full(A));
trials = 100;
addpath('../code')

gaussian = zeros(length(ds),1);
srtt = zeros(length(ds),1);
sse = zeros(length(ds),1);

ds = round(logspace(2,4,7));
for d_idx = 1:length(ds)
    d = ds(d_idx)

    for trial = 1:trials
        S = randn(d,n) / sqrt(d);
        SQ = S*Q;
        svals = svd(SQ);
        gaussian(d_idx) = gaussian(d_idx) + max(max(svals)-1,1-min(svals))/trials;
    
        signs = 2*randi(2,n,1)-3; s = randsample(n,d);
        SQ = dct(signs .* Q); SQ = sqrt(n/d) * SQ(s,:);
        svals = svd(SQ);
        srtt(d_idx) = srtt(d_idx) + max(max(svals)-1,1-min(svals))/trials;
    
        S = sparsesign(d,n,8);
        SQ = S*Q;
        svals = svd(SQ);
        sse(d_idx) = sse(d_idx) + max(max(svals)-1,1-min(svals))/trials;
    end
end

%% Plot
close all
loglog(ds,gaussian,'^:','LineWidth',2,'MarkerFaceColor',"#0072BD",'MarkerSize',16); hold on
loglog(ds,srtt,'s:','LineWidth',2,'MarkerFaceColor',"#D95319",'MarkerSize',16)
loglog(ds,sse,'o:','LineWidth',2,'MarkerFaceColor',"#EDB120",'MarkerSize',16)
loglog(ds,sqrt(k./ds),'k--','LineWidth',3)
xlabel('Sketching dimension $d$')
ylabel('Distortion $\varepsilon$')
legend({'Gaussian', 'SRTT', 'Sparse sign', '$\sqrt{k/d}$'})
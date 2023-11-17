n = 1e6;
k = 200;
d = 2*k;
b = randn(n,1);
A = randn(n,k);
trials = 20;
fprintf('construction / apply vector / apply matrix\n')
addpath('../code')

% Gaussian embedding
construction = 0;
applyvec = 0;
applymat = 0;
tic; 
for i = 1:trials
    tic; S = randn(d,n)/sqrt(d); construction = construction + toc/trials;
    tic; Sb = S*b; applyvec = applyvec + toc/trials;
    tic; SA = S*A; applymat = applymat + toc/trials;
end
fprintf('Gaussian embedding: %f / %f / %f\n',construction,applyvec,applymat)

% SRDCT
construction = 0;
applyvec = 0;
applymat = 0;
tic; 
for i = 1:trials
    tic; signs = 2*randi(2,n,1)-3; s = randsample(n,d); construction = construction + toc/trials;
    tic; Sb = dct(signs .* b); Sb = sqrt(n/d) * Sb(s,:); applyvec = applyvec + toc/trials;
    tic; SA = dct(signs .* A); SA = sqrt(n/d) * SA(s,:); applymat = applymat + toc/trials;
end
fprintf('SRTT: %f / %f / %f\n',construction,applyvec,applymat)

% Gaussian embedding
construction = 0;
applyvec = 0;
applymat = 0;
tic; 
for i = 1:trials
    tic; S = sparsesign(d,n,8); construction = construction + toc/trials;
    tic; Sb = S*b; applyvec = applyvec + toc/trials;
    tic; SA = S*A; applymat = applymat + toc/trials;
end
fprintf('Sparse sign embedding: %f / %f / %f\n',construction,applyvec,applymat)
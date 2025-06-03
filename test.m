N = 1000;
A = sprandsym(N, 0.01, 0.1, 1);

M = 10;
sigma = zeros(M, 1);
for m = 0:M-1
    theta = 2 * pi * (m + 0.5) / M;
    sigma(m+1) = 1.0 + 0.1*exp(1i * theta);
end

b = ones(N, 1);

max_itr = 100000;
threshold = 1e-13;

[x_CG, flag_CG, rres_CG, itrs_CG] = shifted_CG(A, b, N, sigma, M, max_itr, threshold);
[x_MR, flag_MR, rres_MR, itrs_MR] = shifted_MINRES(A, b, N, sigma, M, max_itr, threshold);

% 実行結果の検証
true_res_CG = zeros(M,1);
true_res_MR = zeros(M,1);
for m = 1:1:M
    r = b - (A*x_CG(:,m) + sigma(m)*x_CG(:,m));
    true_res_CG(m) = norm(r)/norm(b);
    r = b - (A*x_MR(:,m) + sigma(m)*x_MR(:,m));
    true_res_MR(m) = norm(r)/norm(b);
end

%
% verify program for my implementations
% First update : 2024/12/13
% Last update  : 2024/12/13
% Created by "ShunHidaka (https://github.com/ShunHidaka)"
%

[A, rows, cols, entries] = mmread("ELSES_MATRIX_VCNT4000std_A.mtx");
N = rows;
sigma = 1.0 + 0.001;
b = ones(N, 1);
max_itr = 100000;
threshold = 1e-13;

[x_cg,flag_cg,relres_cg,iter_cg]  = pcg(A+sigma*eye(N), b, threshold, max_itr);
[x_mr,flag_mr,relres_mr,iter_mr]  = minres(A+sigma*eye(N), b, threshold, max_itr);

[x_my_mr, flag_my_mr, relres_my_mr, iter_my_mr] = shifted_MINRES(A, b, N, sigma, 1, max_itr, threshold);
[x_my_cg, flag_my_cg, relres_my_cg, iter_my_cg] = shifted_CG(A, b, N, sigma, 1, max_itr, threshold);
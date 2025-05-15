%
% sample program for shifted Krylov solver
% First update : 2024/12/12
% Last update  : 2024/12/13
% Created by "ShunHidaka (https://github.com/ShunHidaka)"
%

% 行列Aの用意
% https://math.nist.gov/MatrixMarket/mmio/matlab/mmiomatlab.html
% http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_VCNT4000std_20130517.tgz
%[A, rows, cols, entries] = mmread("ELSES_MATRIX_VCNT4000std_A.mtx");
%[A, rows, cols, entries] = mmread("ELSES_MATRIX_VCNT900_A.mtx");
[A, rows, cols, entries] = mmread("ELSES_MATRIX_VCNT900h_A.mtx");
N = rows;
% シフトsigmaの用意
M = 10;
sigma = zeros(M, 1);
for m = 1:1:M
    %sigma(m) = 0.5 + 0.001i*m;
    %sigma(m) = 0.001*m + 0.01i;
    sigma(m) = 0.0 + 0.01 * exp(1i * 2 * pi * (m+0.5) / M);
end
% 右辺ベクトルbの用意
b = ones(N, 1);
% 収束関連の定数を設定
max_itr = 100000;
threshold = 1e-13;

% 解かせる
[x, flag, rres, itrs] = shifted_MINRES(A, b, N, sigma, M, max_itr, threshold);

% 実行結果の検証
true_res = zeros(M,1);
for m = 1:1:M
    r = b - (A*x(:,m) + sigma(m)*x(:,m));
    true_res(m) = norm(r)/norm(b);
end

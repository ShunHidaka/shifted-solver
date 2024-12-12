%
% sample program for shifted Krylov solver
% First update : 2024/12/12
% Last update  : 2024/12/12
% Created by "ShunHidaka (https://github.com/ShunHidaka)"
%

% 行列Aの用意
% https://math.nist.gov/MatrixMarket/mmio/matlab/mmiomatlab.html
[A, rows, cols, entries] = mmread("ELSES_MATRIX_VCNT4000std_A.mtx");
N = rows;
% シフトsigmaの用意
M = 3;
sigma = zeros(M, 1);
for m = 1:1:M
    % VCNT4000std_Aは正定値でないので最小固有値をシフトさせて正定値にする
    sigma(m) = 1.0 + 0.001*m;
end
% 右辺ベクトルbの用意
b = ones(N, 1);
% 収束関連の定数を設定
max_itr = 100000;
threshold = 1e-13;

% 解かせる
x = shifted_CG(A, b, N, sigma, M, max_itr, threshold);

% 実行結果の検証
rel_nrms = zeros(M,1);
for m = 1:1:M
    r = b - (A + sigma(m)*eye(N))*x(:,m);
    rel_nrms(m) = norm(r)/norm(b);
end

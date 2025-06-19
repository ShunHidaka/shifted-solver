% 行列やシフトの設定を変更して検証できるプログラム
% CG法は正定値でなければ動作しないことに注意
% 複数シフトの場合の検証は標準関数が対応していないので
% 自分で繰り返し解かせるなどして行ってください

% 行列とシフトの設定
N = 10000;                % 行列サイズの設定
A = gallery('lehmer', N); % 行列の設定
sigma = 1.0 + 0.001;      % シフトの設定
b = ones(N, 1);           % 右辺ベクトルの設定
max_itr = 1000000;        % 最大反復回数の設定(行列サイズ以上にしておくべき)
threshold = 1e-13;        % 閾値の設定（現在の残差/初期残差 < threshold）
                          % tolerance は threshold よりも大きくすること

% 標準関数での求解
[x_cg,flag_cg,relres_cg,iter_cg]  = pcg(A+sigma*eye(N), b, threshold, max_itr);
[x_mr,flag_mr,relres_mr,iter_mr]  = minres(A+sigma*eye(N), b, threshold, max_itr);

% 自作関数での求解
[x_my_cg, flag_my_cg, relres_my_cg, iter_my_cg] = shifted_CG(A, b, N, sigma, 1, max_itr, threshold);
[x_my_mr, flag_my_mr, relres_my_mr, iter_my_mr] = shifted_MINRES(A, b, N, sigma, 1, max_itr, threshold);

% 結果の比較
% (標準関数の解 - 自作関数の解)で全ての要素の差が
% tolerance 以下なら true(1) が、以上なら false (0) が返る
tolerance = 1e-12;
all(abs(x_cg - x_my_cg) < tolerance)
all(abs(x_mr - x_my_mr) < tolerance)
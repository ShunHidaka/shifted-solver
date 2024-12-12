%
% shifted CG method by MATLAB
% First update : 2024/12/12
% Last update  : 2024/12/12
% Created by "ShunHidaka (https://github.com/ShunHidaka)"
% 実数向けに作成しているため内積の計算に注意！複素数用じゃない
% 引数：
%   行列A, 右辺ベクトルrhs, 行列サイズN,
% 　シフトsigma, シフト数M,
%   最大反復回数max_itr, 閾値threshold
% 返り値：
%   近似解x(M行N列)
%

function x = shifted_CG(A, rhs, N, sigma, M, max_itr, threshold)
    % 使用する変数の宣言
    x       = zeros(N, M); % 近似解
    r       = zeros(N, M); % 残差
    p       = zeros(N, M); % 補助ベクトル
    alpha   = zeros(2, M); % alpha_{j-1}, alpha_{j}
    beta    = zeros(M, 1);
    pi      = zeros(3, M); % pi_{j-1}, pi_{j}, pi_{j+1}
    is_conv = zeros(M, 1);
    Ap      = zeros(N, 1); % 行列ベクトル積で使う一時ベクトル
    rr      = zeros(2, 1); % 演算数を減らすための一時変数 rr_{j-1}, rr_{j}
    % 変数の初期化
    for m = 1:1:M
        r(:,m)     = rhs; % 初期近似解はゼロベクトルとしている
        alpha(2,m) = 1;
        beta(m)    = 0;
        pi(1,m)    = 1;
        pi(2,m)    = 1;
        is_conv(m) = 0;
    end
    r0nrm = norm(rhs);
    % メインループ
    conv_num = 0;
    s = 2;
    rr(2) = dot(r(:,s),r(:,s));
    for j = 1:1:max_itr
        % seed方程式の更新
        p(:,s) = r(:,s) + beta(s)*p(:,s);
        alpha(1,s) = alpha(2,s);
        Ap(:) = A*p(:,s) + sigma(s)*p(:,s);
        alpha(2,s) = rr(2) / dot(p(:,s), Ap(:));
        x(:,s) = x(:,s) + alpha(2,s)*p(:,s);
        r(:,s) = r(:,s) - alpha(2,s)*Ap(:);
        % add方程式の更新
        for m = 1:1:M
            if is_conv(m) ~= 0 || m == s
                continue;
            end
            pi(3,m)    = (1 + (beta(s)/alpha(1,s))*alpha(2,s) + alpha(2,s)*(sigma(m)-sigma(s)))*pi(2,m) - (beta(s)/alpha(1,s))*alpha(2,s)*pi(1,m);
            alpha(2,m) = (pi(2,m)/pi(3,m))   * alpha(2,s);
            beta(m)    = (pi(1,m)/pi(2,m))^2 * beta(s);
            p(:,m) = r(:,m) + beta(m)*p(:,m);
            x(:,m) = x(:,m) + alpha(2,m)*p(:,m);
            r(:,m) = (1.0/pi(3,m))*r(:,s);
            pi(1,m) = pi(2,m);
            pi(2,m) = pi(3,m);
            if norm(r(:,m))/r0nrm <= threshold
                is_conv(m) = j;
                conv_num = conv_num + 1;
            end
        end
        rr(1) = rr(2);
        rr(2) = dot(r(:,s), r(:,s));
        beta(s) = rr(2)/rr(1);
        if norm(r(:,s))/r0nrm <= threshold && is_conv(s) == 0
            is_conv(s) = j;
            conv_num = conv_num + 1;
        end
        % seed switching

        % determine convergence
        if conv_num == M
            break;
        end
    end
    j

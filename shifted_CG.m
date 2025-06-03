%
% shifted CG method by MATLAB
% First update : 2024/12/12
% Last update  : 2025/06/03
% Created by "ShunHidaka (https://github.com/ShunHidaka)"
% 実数向けに作成、複素数でも動く
% 引数：
%   行列A, 右辺ベクトルrhs, 行列サイズN,
% 　シフトsigma, シフト数M,
%   最大反復回数max_itr, 閾値threshold
% 返り値：
%   近似解x（N行M列）
%   フラグflag
%   相対残差rres（M行）
%   収束したときの反復回数itrs（M行）
%

function [x, flag, rres, itrs] = shifted_CG(A, rhs, N, sigma, M, max_itr, threshold)
    x    = zeros(N, M); % 近似解
    flag = 1;           % 正常に収束したかどうかを示すフラグ、0なら正常に実行
    rres = zeros(M, 1); % 各シフト方程式の相対残差
    itrs = zeros(M, 1); % 収束したときの反復回数
    % 使用する変数の宣言
    r       = zeros(N, M); % 残差
    p       = zeros(N, M); % 補助ベクトル
    alpha   = zeros(2, M); % alpha_{j-1}, alpha_{j}
    beta    = zeros(M, 1); %
    pi      = zeros(3, M); % pi_{j-1}, pi_{j}, pi_{j+1}
    is_conv = zeros(M, 1); % 収束した方程式を管理するベクトル
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
    s = 1;
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
            if is_conv(m) == 1 || m == s
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
            % add方程式の収束判定
            rres(m) = norm(r(:,m))/r0nrm;
            if rres(m) <= threshold
                is_conv(m) = 1;
                itrs(m) = j;
                conv_num = conv_num + 1;
            end
        end
        rr(1) = rr(2);
        rr(2) = dot(r(:,s), r(:,s));
        beta(s) = rr(2)/rr(1);
        % seed方程式の収束判定
        rres(s) = norm(r(:,m))/r0nrm;
        if rres(s) <= threshold && is_conv(s) == 0
            is_conv(s) = 1;
            itrs(s) = j;
            conv_num = conv_num + 1;
        end
        % 全ての方程式が収束したかの判定
        if conv_num == M
            flag = 0;
            break;
        end
        % seed switching
        if is_conv(s) ~= 0
            t = s;
            for m = 1:1:M
                if rres(m) > rres(t) && m ~= s
                    t = m;
                end
            end
            beta(t) = (pi(1,t)/pi(2,t))*(pi(1,t)/pi(2,t))*beta(s);
            for m = 1:1:M
                if m == t, continue; end
                pi(1,t) = pi(1,t) / pi(1,t);
                pi(2,t) = pi(2,t) / pi(2,t);
            end
            %fprintf(1, '# SWITCH [%d] to [%d] in %d : %f %f\n', s, t, j, rres(s), rres(t));
            s = t;
        end
    end
    end
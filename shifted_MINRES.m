%
% shifted MINRES method by MATLAB
% First update : 2024/12/12
% Last update  : 2024/12/13
% Created by "ShunHidaka (https://github.com/ShunHidaka)"
% 引数：
%   行列A, 右辺ベクトルrhs, 行列サイズN,
% 　シフトsigma, シフト数M,
%   最大反復回数max_itr, 閾値threshold
% 返り値：
%   近似解x(N行M列)
%   フラグflag
%   相対残差rres（M行）
%   収束したときの反復回数itrs（M行）
%

function [x, flag, rres, itrs] = shifted_MINRES(A, rhs, N, sigma, M, max_itr, threshold)
% shifted_MINRES  Solves multiple shifted linear systems using the shifted MINRES method.
%
%   This function solves the systems:
%       (A + sigma(m) * I) * x(:,m) = rhs     for m = 1,...,M
%   where A is a real symmetric or complex Hermitian matrix.
%
%   The algorithm is based on a shared Lanczos process, and uses complex
%   Givens rotations (BLAS-style ZROTG) for numerical stability.
%
% INPUTS:
%   A         : (N x N) real-symmetric or Hermitian matrix
%   rhs       : (N x 1) right-hand side vector
%   N         : number of rows/columns in A
%   sigma     : (M x 1) vector of shift parameters
%   M         : number of shifted systems
%   max_itr   : maximum number of iterations
%   threshold : convergence tolerance (relative residual)
%
% OUTPUTS:
%   x     : (N x M) approximate solutions
%   flag  : convergence flag (0 if all systems converged, 1 otherwise)
%   rres  : (M x 1) final relative residuals
%   itrs  : (M x 1) number of iterations to convergence for each system
%

    % Initialize outputs
    x    = zeros(N, M);       % Approximate solutions
    flag = 1;                 % Convergence flag, if flag==0 then success
    rres = zeros(M, 1);       % Relative residuals
    itrs = zeros(M, 1);       % Iteration counts

    % Workspace variables
    q       = zeros(N, 3);    % Lanczos basis vectors: q_{j-1}, q_j, q_{j+1}
    alpha   = zeros(1, 1);    % Lanczos scalar alpha_j
    beta    = zeros(2, 1);    % Lanczos scalar beta_{j-1}, beta_j
    r       = zeros(3, M);    % Tridiagonal entries per shift
    G1      = zeros(2, 2, M); % Givens rotation G_{j-2}
    G2      = zeros(2, 2, M); % Givens rotation G_{j-1}
    G3      = zeros(2, 2, M); % Givens rotation G_j
    p       = zeros(N, 3, M); % Auxiliary vectors: p_{j-2}, p_{j-1}, p_j
    f       = ones(M, 1);     % Scaling factor for each system
    h       = ones(M, 1);     % Residual norms per shift
    is_conv = zeros(M, 1);    % Convergence flag for each system

    % Initialization
    r0nrm = norm(rhs);       % Initial residual norm
    q(:,2) = rhs / r0nrm;    % First Lanczos vector q_1
    conv_num = 0;            % Number of converged systems
    h = r0nrm * h;           % Initiali residual norm in Algorithm

    % Main iteration loop
    for j = 1:1:max_itr
        % === Lanczos step ===
        q(:,3) = A*q(:,2) - beta(1)*q(:,1);
        alpha(1) = real(q(:,2)' * q(:,3));
        q(:,3) = q(:,3) - alpha(1)*q(:,2);
        beta(2) = norm(q(:,3));

        % === Update each shifted system ===
        for m = 1:1:M
            if is_conv(m) == 1
                continue;
            end

            % Tridiagonal matrix entry set
            r(1,m)=0;
            r(2,m)=beta(1);
            r(3,m)=alpha(1)+sigma(m);

            % Apply and Generate Givens rotation
            if j >= 3
                r(1:2,m) = G1(:,:,m)*r(1:2,m);
            end
            if j >= 2
                r(2:3,m) = G2(:,:,m)*r(2:3,m);
            end
            G3(:,:,m) = my_zrotg(r(3,m), beta(2));
            r(3,m) = G3(1,:,m) * [r(3,m); beta(2)];

            % Update auxiliary vector p and solution x
            p(:,1,m) = p(:,2,m);
            p(:,2,m) = p(:,3,m);
            p(:,3,m) = ( q(:,2) - r(1,m)*p(:,1,m) - r(2,m)*p(:,2,m) ) / r(3,m);
            x(:,m) = x(:,m) + r0nrm*G3(1,1,m)*f(m)*p(:,3,m);

            % Update residual estimate
            f(m) = G3(2,1,m)*f(m);
            h(m) = abs(G3(2,1,m))*h(m);

            % Save current Givens for next iteration
            G1(:,:,m) = G2(:,:,m);
            G2(:,:,m) = G3(:,:,m);

            % Check Convergence
            rres(m) = h(m)/r0nrm;
            if rres(m) <= threshold
                is_conv(m) = 1;
                itrs(m) = j;
                conv_num = conv_num + 1;
            end
        end

        % === Update Lanczos basis vectors ===
        q(:,3) = q(:,3)/beta(2);
        q(:,1) = q(:,2);
        q(:,2) = q(:,3);
        beta(1) = beta(2);

        % Exit if all systems have converged
        if conv_num == M
            flag = 0;
            break;
        end
    end
end

%-----------------------------------------------------------
function G = my_zrotg(a, b)
% my_zrotg  Generate complex Givens rotation matrix (BLAS-style)
%
%   G = zrotg(a, b) returns a 2x2 unitary matrix G such that:
%       G * [a; b] = [r; 0]
%   The matrix G is of the form:
%       G = [  c      s
%            -conj(s) c ]
%   where c is real, s is complex, and c^2 + |s|^2 = 1.

    if a == 0 && b == 0
        c = 1;
        s = 0;
    else
        norm_ab = sqrt(abs(a)^2 + abs(b)^2);
        sgn_a = a / abs(a);
        if a == 0
            sgn_a = 1;
        end
        c = abs(a) / norm_ab;
        s = sgn_a * conj(b) / norm_ab;
    end
    G = [c, s; -conj(s), c];
end

function [sol, K, Y, L, diagnostics] = solve_sdp(data, opts)


if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'verbose'),   opts.verbose   = 0;     end
if ~isfield(opts,'solver'),    opts.solver    = 'mosek'; end

% -------------------------
% Unpack data
% -------------------------
B = data.B;
Z     = data.Z;     % n×T
Xm     = data.X;    % n×T
Um     = data.U;    % m×T  (m=1 でも可)

Gamma_Matrix = [Um;Xm];

Phi11  = data.Phi11;
Phi12  = data.Phi12;
Phi22  = data.Phi22;

[n, T] = size(Z);
m = size(Um,1);

% -------------------------
% Decision variables
% -------------------------
% 論文でのXはここではXmと混同しないためにLとする
Y      = sdpvar(n,n,'symmetric');
L      = sdpvar(m,n,'full');
alpha  = sdpvar(1,1);
beta   = sdpvar(1,1);
tDelta = sdpvar(1,1);     % = δ²

% -------------------------
G = [ eye(n)  ,  Z-B*Um ;           % 1: n
    zeros(n,n) , -Xm ;           % 2: n
    zeros(n,n) , zeros(n,T) ;           % 3: m
    zeros(m,n+T)];                 % 7: m

Phi  = [Phi11 Phi12; Phi12' Phi22];

% -------------------------
% LMI blocks (25)(26)(27)
% -------------------------
[F1, F2, F3] = original_thesis.build_lmi_blocks(Y,L,alpha,beta,tDelta,G,Phi,data);

% -------------------------
% Constraints
% -------------------------

tolerance = 1e-8;

% 制約を個別に定義（dual取得のため）
const_1 = [F1 - alpha*F2 >= 0];  % Lambda1に対応 (F2 = G*Phi*G'はalphaを含まない)
const_2 = [F3 >= tolerance*eye(n+m)];  % Lambda3に対応
const_alpha = [alpha >= 0];  % Lambda_alphaに対応
const_beta = [beta >= tolerance];  % Lambda_betaに対応
const_tDelta_lower = [tDelta >= tolerance];  % Lambda_tDelta_lowerに対応
const_Y = [Y >= tolerance*eye(n)];  % Yの制約（必要に応じて）
const_tDelta_upper = [(1 - tolerance)^2 >= tDelta];  % tDeltaの上界制約

constr  = [];
constr  = [constr, const_1];
constr  = [constr, const_2];
constr  = [constr, const_alpha];
constr  = [constr, const_beta];
constr  = [constr, const_tDelta_lower];
constr  = [constr, const_Y];
constr  = [constr, const_tDelta_upper];

% Objective: maximize δ  <=>  minimize tDelta
% pi_matrix = eye(T) - pinv(Gamma_Matrix)*Gamma_Matrix;
obj = -tDelta + 1e-6*norm(Y, 'fro') + 1e-6*norm(L, 'fro');

params = sdpsettings('solver', opts.solver, 'verbose', opts.verbose);
diagnostics = optimize(constr, obj, params);

% -------------------------
% 最適化結果のチェック
% -------------------------
sol.status = diagnostics.problem;

if diagnostics.problem ~= 0
    % 解けなかった場合はエラーを出す
    error('solve_sdp:OptimizationFailed', ...
        '最適化が解けませんでした。status: %d, info: %s', ...
        diagnostics.problem, diagnostics.info);
end

% -------------------------
% Output pack
% -------------------------
sol.Y         = value(Y);
sol.L         = value(L);
sol.alpha     = value(alpha);
sol.beta      = value(beta);
% tDelta = δ^{2}なので、δ = (tDelta)^(1/2)
delta = sqrt(value(tDelta));
sol.delta     = delta;
sol.rho = (1-delta)/(1+delta);  % Theorem 3: ρ = (1-δ)/(1+δ)
sol.objective = value(obj);

K = sol.L/sol.Y;
sol.K = K;

% -------------------------
% Dual values (Lagrange multipliers)
% -------------------------
% KKT条件について：
% - 制約 g(x) >= 0 に対して、ラグランジュ乗数 λ >= 0
% - 相補性条件: λ * g(x) = 0
% - YALMIPのdual関数は、制約 g(x) >= 0 に対して λ >= 0 を返す
%
% 制約の順序に注意：const_1, const_2, const_alpha, const_beta, const_tDelta_lower, const_Y, const_tDelta_upper
try
    sol.Lambda1 = dual(const_1);  % F1 - α*F2 >= 0 のdual（サイズ: (3*n+m) × (3*n+m)）
    % 相補性: Lambda1 * (F1 - α*F2) = 0
    sol.Lambda3 = dual(const_2);  % F3 >= tolerance*eye(n+m) のdual（サイズ: (n+m) × (n+m)）
    % 相補性: Lambda3 * (F3 - tolerance*I) = 0
    sol.Lambda_alpha = dual(const_alpha);  % alpha >= 0 のdual（スカラー）
    % 相補性: Lambda_alpha * alpha = 0
    sol.Lambda_beta = dual(const_beta);  % beta >= tolerance のdual（スカラー）
    % 相補性: Lambda_beta * (beta - tolerance) = 0
    sol.Lambda_tDelta = dual(const_tDelta_lower);  % tDelta >= tolerance のdual（スカラー）
    % 相補性: Lambda_tDelta * (tDelta - tolerance) = 0
    % 注意: const_tDelta_upperのdualは通常は0（上界制約がactiveでない限り）
    sol.Lambda_Y = dual(const_Y);  % Y >= tolerance*eye(n) のdual（サイズ: n × n）
    % 相補性: Lambda_Y * (Y - tolerance*eye(n)) = 0
    sol.Lambda_tDelta_upper = dual(const_tDelta_upper);  % (1 - tolerance)^2 >= tDelta のdual（スカラー）
    % 相補性: Lambda_tDelta_upper * ((1 - tolerance)^2 - tDelta) = 0
catch ME
    % dual取得に失敗した場合はエラーを出す（適当な値を返さない）
    error('solve_sdp:DualFailed', ...
        'Dual変数の取得に失敗しました: %s', ME.message);
end
end

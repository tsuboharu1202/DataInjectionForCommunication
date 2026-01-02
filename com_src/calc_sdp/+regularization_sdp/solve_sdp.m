function [sol, K, delta_val, L_val, diagnostics] = solve_sdp(data, gamma, opts)

if nargin < 2 || isempty(gamma), gamma = 1e5; end
if nargin < 3 || isempty(opts), opts = struct(); end
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


[n, T] = size(Z);
m = size(Um,1);

% -------------------------
% Decision variables
% -------------------------
% 論文でのXはここではXmと混同しないためにLとする
L      = sdpvar(T,n,'full');
delta = sdpvar(1,1);     % = δ²


% -------------------------
% LMI blocks
% -------------------------
Xm_L_Sym = (Xm*L + (Xm*L)')/2;
const_mat = [Xm_L_Sym, (Z*L)', zeros(n,m),(Um*L)';
    Z*L, Xm_L_Sym, delta*B, zeros(n,m);
    zeros(m,n), delta*B', eye(m), zeros(m,m);
    Um*L, zeros(m,n), zeros(m,m), eye(m)];

% -------------------------
% Constraints
% -------------------------

tolerance = 1e-8;

% 制約を個別に定義（dual取得のため）
const_1 = [const_mat >= 0];
const_2 = [Xm_L_Sym >= tolerance*eye(n)];
const_3 = [Xm*L == (Xm*L)'];

constr  = [];
constr  = [constr, const_1];
constr  = [constr, const_2];
constr  = [constr, const_3];
% constr  = [constr, delta >= tolerance];
% Objective: maximize δ  <=>  minimize tDelta
pi_matrix = eye(T) - pinv(Gamma_Matrix)*Gamma_Matrix;
obj = -delta + gamma*norm(pi_matrix*L, 'fro');

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
L_val = value(L);
delta_val = value(delta);
sol.L = L_val;
rho_val = (1-delta_val)/(1+delta_val);
sol.rho = rho_val;
sol.objective = value(obj);
sol.delta = delta_val;
sol.F1 = value(const_mat);

K = Um*L_val*(Xm*L_val)^(-1);
sol.K = K;
% -------------------------
% Dual values (Lagrange multipliers)
% -------------------------
try
    sol.Lambda1 = dual(const_1);
    sol.Lambda2 = dual(const_2);
    sol.Lambda3 = dual(const_3);
catch ME
    % dual取得に失敗した場合はエラーを出す（適当な値を返さない）
    error('solve_sdp:DualFailed', ...
        'Dual変数の取得に失敗しました: %s', ME.message);
end
end

% function [sol, K, delta_val, L_val, diagnostics] = solve_sdp_with_robust(data, epsilon, gamma, opts)
% % solve_sdp_with_robust: ロバスト制約付き正則化SDP
% % epsilon項を含む制約行列を使用
% %
% % 引数:
% %   data: SystemData object
% %   epsilon: ロバスト制約パラメータ
% %   gamma: 正則化パラメータ（デフォルト: 1e5）
% %   opts: オプション

% if nargin < 3 || isempty(gamma), gamma = 1e5; end
% if nargin < 4 || isempty(opts), opts = struct(); end
% if ~isfield(opts,'verbose'),   opts.verbose   = 0;     end
% if ~isfield(opts,'solver'),    opts.solver    = 'mosek'; end

% % -------------------------
% % Unpack data
% % -------------------------
% B = data.B;
% Z     = data.Z;     % n×T
% Xm     = data.X;    % n×T
% Um     = data.U;    % m×T  (m=1 でも可)

% Gamma_Matrix = [Um;Xm];


% [n, T] = size(Z);
% m = size(Um,1);

% % -------------------------
% % Decision variables
% % -------------------------
% % 論文でのXはここではXmと混同しないためにLとする
% L      = sdpvar(T,n,'full');
% delta = sdpvar(1,1);     % = δ²


% % -------------------------
% % LMI blocks (epsilon項あり)
% % -------------------------
% Xm_L_Sym = (Xm*L + (Xm*L)')/2;

% const_mat = [Xm_L_Sym, (Z*L)', (Um*L)',epsilon*(Xm*L)';
%     Z*L, Xm_L_Sym-delta*(B*B'), zeros(n,m), zeros(n,n);
%     Um*L, zeros(m,n), eye(m), zeros(m,n);
%     epsilon*(Xm*L), zeros(n,n), zeros(n,m), eye(n)];


% % -------------------------
% % Constraints
% % -------------------------

% tolerance = 1e-6;

% % 制約を個別に定義（dual取得のため）
% const_1 = [const_mat >= tolerance*eye(3*n+m)];
% const_2 = [Xm_L_Sym >= (tolerance)*eye(n)];
% const_3 = [Xm*L == (Xm*L)'];

% constr  = [];
% constr  = [constr, const_1];
% constr  = [constr, const_2];
% constr  = [constr, const_3];

% % Objective: maximize δ  <=>  minimize tDelta
% pi_matrix = eye(T) - pinv(Gamma_Matrix)*Gamma_Matrix;
% obj = -delta + gamma*norm(pi_matrix*L, 'fro')^2;

% params = sdpsettings('solver', opts.solver, 'verbose', opts.verbose);
% diagnostics = optimize(constr, obj, params);

% % -------------------------
% % 最適化結果のチェック
% % -------------------------
% sol.status = diagnostics.problem;

% if diagnostics.problem ~= 0
%     % 解けなかった場合はエラーを出す
%     error('solve_sdp_with_robust:OptimizationFailed', ...
%         '最適化が解けませんでした。status: %d, info: %s', ...
%         diagnostics.problem, diagnostics.info);
% end

% % -------------------------
% % Output pack
% % -------------------------
% L_val = value(L);
% delta_val = sqrt(value(delta));
% sol.L = L_val;
% rho_val = (1-delta_val)/(1+delta_val);
% sol.rho = rho_val;
% sol.objective = value(obj);
% sol.delta = delta_val;
% sol.F1 = value(const_mat);

% K = Um*L_val*(Xm*L_val)^(-1);
% sol.K = K;




% % -------------------------
% % Dual values (Lagrange multipliers)
% % -------------------------
% try
%     sol.Lambda1 = dual(const_1);
%     sol.Lambda2 = dual(const_2);
%     sol.Lambda3 = dual(const_3);
% catch ME
%     % dual取得に失敗した場合はエラーを出す（適当な値を返さない）
%     error('solve_sdp_with_robust:DualFailed', ...
%         'Dual変数の取得に失敗しました: %s', ME.message);
% end

% end


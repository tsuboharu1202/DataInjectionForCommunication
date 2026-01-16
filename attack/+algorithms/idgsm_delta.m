function [X_adv, Z_adv, U_adv, history] = idgsm_delta(ori_sd, opts)
% idgsm_delta: Iterative Direct Gradient Sign Method for delta attack
%
% 引数:
%   ori_sd: SystemData オブジェクト（攻撃対象のデータ）
%   opts: オプション構造体（全て必須）
%
% 必須オプション:
%   opts.gamma: 正則化パラメータ（SDP計算用）
%   opts.epsilon: ノルム制約の上限
%   opts.alpha: ステップサイズ
%   opts.max_iteration: 最大反復回数
%   opts.direction: 'positive' or 'negative'
%
% オプション:
%   opts.save_history: 履歴を保存するか（デフォルト: false）
%   opts.save_grad_norms: 勾配ノルムを保存するか（デフォルト: false）
%   opts.normalize_grad: 勾配を正規化するか（デフォルト: true）
%
% 出力:
%   X_adv, Z_adv, U_adv: 攻撃後のデータ
%   history: 履歴（opts.save_history=true の場合のみ）

% === 引数チェック ===
if nargin < 2
    error('idgsm_delta:MissingOpts', 'opts は必須です。');
end

% 必須パラメータのチェック
required_fields = {'gamma', 'epsilon', 'alpha', 'max_iteration', 'direction'};
for i = 1:length(required_fields)
    if ~isfield(opts, required_fields{i})
        error('idgsm_delta:MissingField', ...
            'opts.%s は必須です。', required_fields{i});
    end
end

gamma = opts.gamma;
epsilon = opts.epsilon;
alpha = opts.alpha;
max_iteration = opts.max_iteration;
direction = opts.direction;

% オプションのデフォルト値（履歴・正規化関連のみ）
if ~isfield(opts, 'save_history'),    opts.save_history = false; end
if ~isfield(opts, 'save_grad_norms'), opts.save_grad_norms = false; end
if ~isfield(opts, 'normalize_grad'),  opts.normalize_grad = true; end

save_history = opts.save_history;
save_grad_norms = opts.save_grad_norms;
normalize_grad = opts.normalize_grad;

% === 初期化 ===
iter = 0;
is_continue = true;

dX = zeros(size(ori_sd.X));
dZ = zeros(size(ori_sd.Z));
dU = zeros(size(ori_sd.U));
current_sd = ori_sd;

% 方向の決定
if strcmpi(direction, 'positive')
    dir_mult = 1;   % deltaを大きくする方向
else
    dir_mult = -1;  % deltaを小さくする方向
end

% 履歴の初期化
if save_history
    history = struct();
    history.dX_history = {};
    history.dZ_history = {};
    history.dU_history = {};
    history.delta_history = [];
    history.iter_count = 0;
    
    if save_grad_norms
        history.grad_norms_X_Linf = [];
        history.grad_norms_X_Fro = [];
        history.grad_norms_Z_Linf = [];
        history.grad_norms_Z_Fro = [];
        history.grad_norms_U_Linf = [];
        history.grad_norms_U_Fro = [];
        history.grad_max_ratio_X = [];
        history.grad_max_ratio_Z = [];
        history.grad_max_ratio_U = [];
        history.grad_top3_ratio_X = [];
        history.grad_top3_ratio_Z = [];
        history.grad_top3_ratio_U = [];
    end
else
    history = [];
end

% === メインループ ===
while is_continue
    fprintf('iter: %d\n', iter);
    
    % 勾配を計算
    grad_opts = struct('gamma', gamma);
    [X_grad, Z_grad, U_grad] = common.calc_grad(current_sd, grad_opts);
    
    % 勾配ノルムを保存
    if save_history && save_grad_norms
        history = save_gradient_norms(history, X_grad, Z_grad, U_grad);
    end
    
    % 正規化
    if normalize_grad
        [X_grad, Z_grad, U_grad] = normalize_gradients(X_grad, Z_grad, U_grad);
    end
    
    % ノイズを累積
    dX = dX + alpha * X_grad * dir_mult;
    dZ = dZ + alpha * Z_grad * dir_mult;
    dU = dU + alpha * U_grad * dir_mult;
    
    % 投影
    checkX_temp = false(size(ori_sd.X));
    checkZ_temp = false(size(ori_sd.Z));
    checkU_temp = false(size(ori_sd.U));
    [dX, dZ, dU, ~, ~, ~] = algorithms.projector(dX, dZ, dU, checkX_temp, checkZ_temp, checkU_temp, epsilon);
    
    X_adv = ori_sd.X + dX;
    Z_adv = ori_sd.Z + dZ;
    U_adv = ori_sd.U + dU;
    
    % 攻撃後のデータでSDPを解く
    current_sd = datasim.SystemData(ori_sd.A, ori_sd.B, X_adv, Z_adv, U_adv, ...
        ori_sd.Phi11, ori_sd.Phi12, ori_sd.Phi22);
    
    [sol_temp, ~, ~, ~, diagnostics_temp] = proposed.solve_sdp(current_sd, gamma);
    
    if diagnostics_temp.problem ~= 0
        error('idgsm_delta:SDPInfeasible', ...
            'SDPがinfeasibleまたはエラーです (iter=%d, status=%d)', iter, diagnostics_temp.problem);
    end
    
    delta_temp = sol_temp.delta;
    fprintf('delta_temp: %.6e\n', delta_temp);
    
    % 履歴を保存
    if save_history
        history.dX_history{end+1} = dX;
        history.dZ_history{end+1} = dZ;
        history.dU_history{end+1} = dU;
        history.delta_history(end+1) = delta_temp;
    end
    
    % 終了条件
    is_continue = (iter < max_iteration);
    iter = iter + 1;
end

if save_history
    history.iter_count = iter;
end

end

%% ローカル関数

function [X_grad, Z_grad, U_grad] = normalize_gradients(X_grad, Z_grad, U_grad)
    % 要素単位の最大値で正規化
    max_X = max(abs(X_grad(:)));
    max_Z = max(abs(Z_grad(:)));
    max_U = max(abs(U_grad(:)));
    
    if max_X > eps
        X_grad = X_grad / max_X;
    end
    if max_Z > eps
        Z_grad = Z_grad / max_Z;
    end
    if max_U > eps
        U_grad = U_grad / max_U;
    end
end

function history = save_gradient_norms(history, X_grad, Z_grad, U_grad)
    % L∞ノルム
    max_X = max(abs(X_grad(:)));
    max_Z = max(abs(Z_grad(:)));
    max_U = max(abs(U_grad(:)));
    
    % Frobeniusノルム
    fro_X = norm(X_grad, 'fro');
    fro_Z = norm(Z_grad, 'fro');
    fro_U = norm(U_grad, 'fro');
    
    history.grad_norms_X_Linf(end+1) = max_X;
    history.grad_norms_X_Fro(end+1) = fro_X;
    history.grad_norms_Z_Linf(end+1) = max_Z;
    history.grad_norms_Z_Fro(end+1) = fro_Z;
    history.grad_norms_U_Linf(end+1) = max_U;
    history.grad_norms_U_Fro(end+1) = fro_U;
    
    % 要素分布
    abs_X = abs(X_grad(:));
    abs_Z = abs(Z_grad(:));
    abs_U = abs(U_grad(:));
    
    sum_X = sum(abs_X);
    sum_Z = sum(abs_Z);
    sum_U = sum(abs_U);
    
    history.grad_max_ratio_X(end+1) = safe_divide(max_X, sum_X);
    history.grad_max_ratio_Z(end+1) = safe_divide(max_Z, sum_Z);
    history.grad_max_ratio_U(end+1) = safe_divide(max_U, sum_U);
    
    sorted_X = sort(abs_X, 'descend');
    sorted_Z = sort(abs_Z, 'descend');
    sorted_U = sort(abs_U, 'descend');
    
    history.grad_top3_ratio_X(end+1) = safe_divide(sum(sorted_X(1:min(3, length(sorted_X)))), sum_X);
    history.grad_top3_ratio_Z(end+1) = safe_divide(sum(sorted_Z(1:min(3, length(sorted_Z)))), sum_Z);
    history.grad_top3_ratio_U(end+1) = safe_divide(sum(sorted_U(1:min(3, length(sorted_U)))), sum_U);
end

function result = safe_divide(numerator, denominator)
    if denominator > eps
        result = numerator / denominator;
    else
        result = 0;
    end
end

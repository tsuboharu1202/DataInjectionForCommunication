function [X_adv, Z_adv, U_adv, history] = idgsm_delta(ori_sd, save_history, alpha, escape_local_min, epsilon, direction, save_grad_norms, opts)
% idgsm_delta: Iterative Direct Gradient Sign Method for delta attack
%
% Inputs:
%   ori_sd: Original system data
%   save_history: (optional) if true, save noise history at each step
%   alpha: (optional) step size for IDGSM. If not provided, uses cfg.Const.IDGSM_ALPHA
%   escape_local_min: (optional) if true, add random noise when stuck in local minimum. If not provided, uses cfg.Const.IDGSM_ESCAPE_LOCAL_MIN
%   epsilon: (optional) norm constraint upper limit. If not provided, uses cfg.Const.ATTACKER_UPPERLIMIT
%   direction: (optional) 'positive' (deltaを大きくする) または 'negative' (deltaを小さくする)
%              デフォルト: 'negative'
%   save_grad_norms: (optional) if true and save_history is true, save gradient norms (L∞ and Frobenius) at each step
%                    デフォルト: false
%   opts: (optional) 追加オプション（gammaなど）
%
% Outputs:
%   X_adv, Z_adv, U_adv: Adversarial data
%   history: (optional) struct with fields:
%       - dX_history: cell array of dX at each step
%       - dZ_history: cell array of dZ at each step
%       - dU_history: cell array of dU at each step
%       - delta_history: array of delta at each step
%       - iter_count: number of iterations
%       - grad_norms_X_Linf: (if save_grad_norms=true) array of L∞ norms of X_grad before normalization
%       - grad_norms_X_Fro: (if save_grad_norms=true) array of Frobenius norms of X_grad before normalization
%       - grad_norms_Z_Linf: (if save_grad_norms=true) array of L∞ norms of Z_grad before normalization
%       - grad_norms_Z_Fro: (if save_grad_norms=true) array of Frobenius norms of Z_grad before normalization
%       - grad_norms_U_Linf: (if save_grad_norms=true) array of L∞ norms of U_grad before normalization
%       - grad_norms_U_Fro: (if save_grad_norms=true) array of Frobenius norms of U_grad before normalization

if nargin < 2
    save_history = false;
end
if nargin < 3 || isempty(alpha)
    alpha = cfg.Const.IDGSM_ALPHA;
end
if nargin < 4 || isempty(escape_local_min)
    % 以前は局所最適解回避のためのランダムノイズ付加に使用していたが、
    % 現在はランダムノイズを加えない方針のため、フラグは保持のみで未使用。
    % escape_local_min = false;
end
if nargin < 5 || isempty(epsilon)
    epsilon = cfg.Const.ATTACKER_UPPERLIMIT;
end
if nargin < 6 || isempty(direction)
    direction = 'negative';
end
if nargin < 7 || isempty(save_grad_norms)
    save_grad_norms = false;
end
if nargin < 8
    opts = struct();
end
% optsからも取得可能にする（後方互換性のため）
if isfield(opts, 'save_grad_norms')
    save_grad_norms = opts.save_grad_norms;
end

iter = 0;

is_continue = true;

dX = zeros(size(ori_sd.X));
dZ = zeros(size(ori_sd.Z));
dU = zeros(size(ori_sd.U));
current_sd = ori_sd;

% Initialize history if requested
if save_history
    history = struct();
    history.dX_history = {};
    history.dZ_history = {};
    history.dU_history = {};
    history.delta_history = [];
    history.iter_count = 0;
    
    % Initialize gradient norm arrays if requested
    if save_grad_norms
        history.grad_norms_X_Linf = [];
        history.grad_norms_X_Fro = [];
        history.grad_norms_Z_Linf = [];
        history.grad_norms_Z_Fro = [];
        history.grad_norms_U_Linf = [];
        history.grad_norms_U_Fro = [];
        % 勾配の要素分布
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

while is_continue
    fprintf('iter: %d\n', iter);
    
    % 勾配を計算（regularizationのみ）
    [X_grad, Z_grad, U_grad] = attack.calc_grad(current_sd, opts);
    
    % 方向を適用
    if strcmpi(direction, 'positive')
        dir_mult = 1;  % deltaを大きくする方向
    else
        dir_mult = -1;  % deltaを小さくする方向
    end
    
    % 正規化の有効/無効をoptsから取得（デフォルトはtrue）
    normalize_grad = true;
    if isfield(opts, 'normalize_grad')
        normalize_grad = opts.normalize_grad;
    end
    
    % 正規化前の勾配ノルムを保存（save_grad_normsがtrueの場合）
    if save_history && save_grad_norms
        % L∞ノルム（要素単位の最大値）
        max_X = max(abs(X_grad(:)));
        max_Z = max(abs(Z_grad(:)));
        max_U = max(abs(U_grad(:)));
        
        % Frobeniusノルム
        fro_X = norm(X_grad, 'fro');
        fro_Z = norm(Z_grad, 'fro');
        fro_U = norm(U_grad, 'fro');
        
        % 保存
        history.grad_norms_X_Linf(end+1) = max_X;
        history.grad_norms_X_Fro(end+1) = fro_X;
        history.grad_norms_Z_Linf(end+1) = max_Z;
        history.grad_norms_Z_Fro(end+1) = fro_Z;
        history.grad_norms_U_Linf(end+1) = max_U;
        history.grad_norms_U_Fro(end+1) = fro_U;
        
        % 勾配の要素分布を確認（最大要素の割合など）
        % 最大要素の絶対値が全体の何%を占めるか
        abs_X = abs(X_grad(:));
        abs_Z = abs(Z_grad(:));
        abs_U = abs(U_grad(:));
        
        sum_X = sum(abs_X);
        sum_Z = sum(abs_Z);
        sum_U = sum(abs_U);
        
        if sum_X > eps
            history.grad_max_ratio_X(end+1) = max_X / sum_X;
        else
            history.grad_max_ratio_X(end+1) = 0;
        end
        if sum_Z > eps
            history.grad_max_ratio_Z(end+1) = max_Z / sum_Z;
        else
            history.grad_max_ratio_Z(end+1) = 0;
        end
        if sum_U > eps
            history.grad_max_ratio_U(end+1) = max_U / sum_U;
        else
            history.grad_max_ratio_U(end+1) = 0;
        end
        
        % 上位3要素の合計が全体の何%を占めるか
        sorted_X = sort(abs_X, 'descend');
        sorted_Z = sort(abs_Z, 'descend');
        sorted_U = sort(abs_U, 'descend');
        
        if sum_X > eps
            history.grad_top3_ratio_X(end+1) = sum(sorted_X(1:min(3, length(sorted_X)))) / sum_X;
        else
            history.grad_top3_ratio_X(end+1) = 0;
        end
        if sum_Z > eps
            history.grad_top3_ratio_Z(end+1) = sum(sorted_Z(1:min(3, length(sorted_Z)))) / sum_Z;
        else
            history.grad_top3_ratio_Z(end+1) = 0;
        end
        if sum_U > eps
            history.grad_top3_ratio_U(end+1) = sum(sorted_U(1:min(3, length(sorted_U)))) / sum_U;
        else
            history.grad_top3_ratio_U(end+1) = 0;
        end
    end
    
    % 正規化（要素単位の最大値で正規化）
    if normalize_grad
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
    
    % ノイズを累積
    dX = dX + alpha * X_grad * dir_mult;
    dZ = dZ + alpha * Z_grad * dir_mult;
    dU = dU + alpha * U_grad * dir_mult;
    
    % 投影（IDGSMでは毎回新しいフラグを使用）
    checkX_temp = false(size(ori_sd.X));
    checkZ_temp = false(size(ori_sd.Z));
    checkU_temp = false(size(ori_sd.U));
    [dX, dZ, dU, ~, ~, ~] = attack.projector(dX, dZ, dU, checkX_temp, checkZ_temp, checkU_temp, epsilon);
    
    X_adv = ori_sd.X + dX;
    Z_adv = ori_sd.Z + dZ;
    U_adv = ori_sd.U + dU;
    
    % 攻撃後のデータでSDPを解く（regularizationのみ）
    current_sd = datasim.SystemData(ori_sd.A, ori_sd.B, X_adv, Z_adv, U_adv, ...
        ori_sd.Phi11, ori_sd.Phi12, ori_sd.Phi22);
    
    gamma = 1e3;
    if isfield(opts, 'gamma')
        gamma = opts.gamma;
    end
    [sol_temp, ~, ~, ~, diagnostics_temp] = regularization_sdp.solve_sdp(current_sd, gamma);
    
    % SDPがinfeasibleまたはエラーの場合はエラーを投げる
    if diagnostics_temp.problem ~= 0
        error('idgsm_delta:SDPInfeasible', ...
            'SDPがinfeasibleまたはエラーです (iter=%d, status=%d)', iter, diagnostics_temp.problem);
    end
    
    delta_temp = sol_temp.delta;
    fprintf('delta_temp: %.6e\n', delta_temp);
    
    % Save history if requested
    if save_history
        history.dX_history{end+1} = dX;
        history.dZ_history{end+1} = dZ;
        history.dU_history{end+1} = dU;
        history.delta_history(end+1) = delta_temp;
    end
    
    % 終了条件: 最大反復回数に達した
    is_continue = (iter < cfg.Const.MAX_ITERATION);
    
    iter = iter + 1;  % Increment the iteration counter
end

if save_history
    history.iter_count = iter;
end

end


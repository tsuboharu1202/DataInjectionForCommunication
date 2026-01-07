function [X_adv, Z_adv, U_adv, history] = idgsm_delta(ori_sd, save_history, alpha, escape_local_min, epsilon, direction, opts)
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

if nargin < 2
    save_history = false;
end
if nargin < 3 || isempty(alpha)
    alpha = cfg.Const.IDGSM_ALPHA;
end
if nargin < 4 || isempty(escape_local_min)
    % 以前は局所最適解回避のためのランダムノイズ付加に使用していたが、
    % 現在はランダムノイズを加えない方針のため、フラグは保持のみで未使用。
    escape_local_min = false;
end
if nargin < 5 || isempty(epsilon)
    epsilon = cfg.Const.ATTACKER_UPPERLIMIT;
end
if nargin < 6 || isempty(direction)
    direction = 'negative';
end
if nargin < 7
    opts = struct();
end

iter = 0;

checkX = false(size(ori_sd.X));   % dX と同じ大きさの logical false
checkZ = false(size(ori_sd.Z));
checkU = false(size(ori_sd.U));

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
    
    % X_grad, Z_grad, U_gradの正規化（要素単位の最大値で正規化）
    % これにより、各勾配の最大要素が1になり、alphaの効果が明確になる
    % また、projectorの要素単位制約（L∞ノルム）と整合性がある
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
    
    % ノイズを累積
    dX = dX + alpha * X_grad * dir_mult;
    dZ = dZ + alpha * Z_grad * dir_mult;
    dU = dU + alpha * U_grad * dir_mult;
    
    % 投影
    [dX, dZ, dU, checkX, checkZ, checkU] = attack.projector(dX, dZ, dU, checkX, checkZ, checkU, epsilon);
    
    okX = all(checkX(:));   % R2018b+ なら all(checkX,"all") でもOK
    okZ = all(checkZ(:));
    okU = all(checkU(:));
    
    X_adv = ori_sd.X + dX;
    Z_adv = ori_sd.Z + dZ;
    U_adv = ori_sd.U + dU;
    
    allDone = okX && okZ && okU;
    
    % 攻撃後のデータでSDPを解く（regularizationのみ）
    current_sd = datasim.SystemData(ori_sd.A, ori_sd.B, X_adv, Z_adv, U_adv, ...
        ori_sd.Phi11, ori_sd.Phi12, ori_sd.Phi22);
    
    gamma = 1e3;
    if isfield(opts, 'gamma')
        gamma = opts.gamma;
    end
    [sol_temp, ~, ~, ~, ~] = regularization_sdp.solve_sdp(current_sd, gamma);
    delta_temp = sol_temp.delta;
    fprintf('delta_temp: %.6e\n', delta_temp);
    
    % Save history if requested
    if save_history
        history.dX_history{end+1} = dX;
        history.dZ_history{end+1} = dZ;
        history.dU_history{end+1} = dU;
        history.delta_history(end+1) = delta_temp;
    end
    
    % 終了条件: すべての要素が処理済み または 最大反復回数に達した または deltaが目標範囲外
    is_continue = ~allDone && (iter < cfg.Const.MAX_ITERATION);
    
    iter = iter + 1;  % Increment the iteration counter
end

if save_history
    history.iter_count = iter;
end

end


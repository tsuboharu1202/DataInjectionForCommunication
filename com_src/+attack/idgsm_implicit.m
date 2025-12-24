function [X_adv, Z_adv, U_adv, history] = idgsm_implicit(sd, eps_att, max_iter, save_history, minimize_delta)
% idgsm_implicit: Iterative Direct Gradient Sign Method using implicit differentiation
%                 to maximize or minimize Delta (rho)
%
% This function uses implicit differentiation to compute gradients of Delta,
% which is related to the spectral radius rho through a monotonic relationship.
%
% Inputs:
%   sd: SystemData object
%   eps_att: (optional) Attack parameter (default: from cfg.Const.ATTACKER_UPPERLIMIT)
%   max_iter: (optional) Maximum number of iterations (default: 10)
%   save_history: (optional) If true, save attack history (default: false)
%   minimize_delta: (optional) If true, minimize Delta (use -dtDelta/dD).
%                    If false (default), maximize Delta (use dtDelta/dD).
%
% Outputs:
%   X_adv, Z_adv, U_adv: Adversarial data
%   history: (optional) Attack history containing rho values at each iteration

if nargin < 2 || isempty(eps_att)
    eps_att = cfg.Const.ATTACKER_UPPERLIMIT;
end
if nargin < 3 || isempty(max_iter)
    max_iter = 10;
end
if nargin < 4
    save_history = false;
end
if nargin < 5
    minimize_delta = false;
end

% 初期データ
X_adv = sd.X;
Z_adv = sd.Z;
U_adv = sd.U;

% 履歴を保存する場合
if save_history
    history = struct();
    history.rho = zeros(max_iter + 1, 1);
    
    % 初期rhoを計算
    sd_init = datasim.SystemData(sd.A, sd.B, X_adv, Z_adv, U_adv, sd.Phi11, sd.Phi12, sd.Phi22);
    [sol_init, K_init, ~, ~, ~] = original_thesis.solve_sdp(sd_init);
    Ac_init = sd.A + sd.B * K_init;
    lambda_init = eig(Ac_init);
    history.rho(1) = max(abs(lambda_init));
end

% 反復的に攻撃を実行
for iter = 1:max_iter
    % 現在のデータでSystemDataを作成
    sd_current = datasim.SystemData(sd.A, sd.B, X_adv, Z_adv, U_adv, sd.Phi11, sd.Phi12, sd.Phi22);
    
    % Implicit differentiationで勾配を計算
    try
        [gradX, gradZ, gradU] = attack.calc_grad_implicit(sd_current, minimize_delta);
    catch ME
        warning('attack:idgsm:gradError', '勾配計算でエラーが発生しました（反復 %d）: %s', iter, ME.message);
        break;
    end
    
    % 勾配の符号に基づいてデータを更新
    % d_adv = d_ori + eps * |d_ori| * sign(real(grad))
    if ~isempty(X_adv)
        X_adv = X_adv + eps_att * abs(X_adv) .* sign(real(gradX));
    end
    if ~isempty(Z_adv)
        Z_adv = Z_adv + eps_att * abs(Z_adv) .* sign(real(gradZ));
    end
    if ~isempty(U_adv)
        U_adv = U_adv + eps_att * abs(U_adv) .* sign(real(gradU));
    end
    
    % 履歴を保存する場合
    if save_history
        sd_updated = datasim.SystemData(sd.A, sd.B, X_adv, Z_adv, U_adv, sd.Phi11, sd.Phi12, sd.Phi22);
        try
            [sol_updated, K_updated, ~, ~, ~] = original_thesis.solve_sdp(sd_updated);
            Ac_updated = sd.A + sd.B * K_updated;
            lambda_updated = eig(Ac_updated);
            history.rho(iter + 1) = max(abs(lambda_updated));
        catch
            history.rho(iter + 1) = NaN;
        end
    end
end

% 履歴を保存しない場合
if ~save_history
    history = [];
end

end


function [X_adv, Z_adv, U_adv, history] = execute_attack(system_data, method, eps_att, direction, save_history, opts)
% execute_attack: Execute attack method on system data
%
% Inputs:
%   system_data: System data to attack
%   method: Attack method (cfg.AttackType)
%   eps_att: (optional) 攻撃強度（デフォルト: cfg.Const.ATTACKER_UPPERLIMIT）
%   direction: (optional) 'positive' (deltaを大きくする) または 'negative' (deltaを小さくする)
%              デフォルト: 'negative'
%   save_history: (optional) If true, save attack history (only for IDGSM methods)
%   opts: (optional) 追加オプション（gamma、alpha、escape_local_minなど）
%
% Outputs:
%   X_adv, Z_adv, U_adv: Adversarial data
%   history: (optional) Attack history（IDGSM methodsの場合のみ）

if nargin < 3
    eps_att = [];
end
if nargin < 4 || isempty(direction)
    direction = 'negative';
end
if nargin < 5
    save_history = false;
end
if nargin < 6
    opts = struct();
end

history = [];

% methodを文字列に変換
method_str = string(method);

switch method_str
    case {cfg.AttackType.DIRECT_DGSM_DELTA, "DIRECT_DGSM_DELTA"}
        [X_adv, Z_adv, U_adv] = attack.dgsm_delta(system_data, eps_att, direction, opts);
        
    case {cfg.AttackType.IMPLICIT_IDGSM_DELTA, "IMPLICIT_IDGSM_DELTA"}
        % save_grad_normsをoptsから取得（デフォルトはfalse）
        save_grad_norms = false;
        if isfield(opts, 'save_grad_norms')
            save_grad_norms = opts.save_grad_norms;
        end
        if save_history
            [X_adv, Z_adv, U_adv, history] = attack.idgsm_delta(system_data, save_history, [], [], eps_att, direction, save_grad_norms, opts);
        else
            [X_adv, Z_adv, U_adv] = attack.idgsm_delta(system_data, false, [], [], eps_att, direction, save_grad_norms, opts);
        end
        
    case {cfg.AttackType.IMPLICIT_IDGSM_DELTA_POSITIVE, "IMPLICIT_IDGSM_DELTA_POSITIVE"}
        % deltaを大きくする方向
        % save_grad_normsをoptsから取得（デフォルトはfalse）
        save_grad_norms = false;
        if isfield(opts, 'save_grad_norms')
            save_grad_norms = opts.save_grad_norms;
        end
        if save_history
            [X_adv, Z_adv, U_adv, history] = attack.idgsm_delta(system_data, save_history, [], [], eps_att, 'positive', save_grad_norms, opts);
        else
            [X_adv, Z_adv, U_adv] = attack.idgsm_delta(system_data, false, [], [], eps_att, 'positive', save_grad_norms, opts);
        end
        
    case {cfg.AttackType.IMPLICIT_IDGSM_DELTA_NEGATIVE, "IMPLICIT_IDGSM_DELTA_NEGATIVE"}
        % deltaを小さくする方向
        % save_grad_normsをoptsから取得（デフォルトはfalse）
        save_grad_norms = false;
        if isfield(opts, 'save_grad_norms')
            save_grad_norms = opts.save_grad_norms;
        end
        if save_history
            [X_adv, Z_adv, U_adv, history] = attack.idgsm_delta(system_data, save_history, [], [], eps_att, 'negative', save_grad_norms, opts);
        else
            [X_adv, Z_adv, U_adv] = attack.idgsm_delta(system_data, false, [], [], eps_att, 'negative', save_grad_norms, opts);
        end
        
    otherwise
        error('attack:unknownMethod','Unknown attack method: %s', method_str);
end
end
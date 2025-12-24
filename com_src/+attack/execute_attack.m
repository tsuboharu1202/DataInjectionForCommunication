function [X_adv, Z_adv, U_adv, history] = execute_attack(system_data, method, eps_att, max_iter, save_history)
% execute_attack: Execute attack method on system data
%
% Inputs:
%   system_data: System data to attack
%   method: Attack method (cfg.AttackType)
%   eps_att: (optional) Attack parameter (for DGSM methods)
%   max_iter: (optional) Maximum number of iterations (only for IDGSM methods)
%   save_history: (optional) If true, save attack history (only for IDGSM methods)
%
% Outputs:
%   X_adv, Z_adv, U_adv: Adversarial data
%   history: (optional) Attack history (only for IDGSM methods if save_history=true)

if nargin < 3
    eps_att = [];
end
if nargin < 4
    max_iter = [];
end
if nargin < 5
    save_history = false;
end

history = [];

switch method
    case cfg.AttackType.DIRECT_DGSM_EV
        [X_adv, Z_adv, U_adv] = attack.dgsm_ev(system_data, true, eps_att);
        
    case cfg.AttackType.IMPLICIT_IDGSM_EV
        [X_adv, Z_adv, U_adv, history] = attack.idgsm_implicit(system_data, eps_att, max_iter, save_history, false);
        
    case cfg.AttackType.IMPLICIT_IDGSM_RHO_LARGE
        % 勾配と同じ方向にノイズを加える（Deltaを大きくする、rhoを大きくする）
        [X_adv, Z_adv, U_adv, history] = attack.idgsm_rho_large(system_data, eps_att, max_iter, save_history);
        
    case cfg.AttackType.IMPLICIT_IDGSM_RHO_SMALL
        % 勾配と逆方向にノイズを加える（Deltaを小さくする、rhoを小さくする）
        [X_adv, Z_adv, U_adv, history] = attack.idgsm_rho_small(system_data, eps_att, max_iter, save_history);
        
    otherwise
        error('attack:unknownMethod','Unknown attack method: %s', string(method));
end
end
function [X_adv, Z_adv, U_adv, history] = execute_attack(system_data, method, eps_att, save_history)
% execute_attack: Execute attack method on system data
%
% Inputs:
%   system_data: System data to attack
%   method: Attack method (cfg.AttackType)
%   eps_att: (optional) Attack parameter (for DGSM methods)
%   save_history: (optional) If true, save attack history (only for IDGSM methods)
%
% Outputs:
%   X_adv, Z_adv, U_adv: Adversarial data
%   history: (optional) Attack history (only for IDGSM methods if save_history=true)

if nargin < 3
    eps_att = [];
end
if nargin < 4
    save_history = false;
end

history = [];

switch method
    case cfg.AttackType.DIRECT_DGSM_EV
        [X_adv, Z_adv, U_adv] = attack.dgsm_ev(system_data, true, eps_att);
        
        % case cfg.AttackType.Direct_IDGSM_EV
        %     % For now, IMPLICIT_IDGSM_EV uses the same implementation as DIRECT_DGSM_EV
        %     % TODO: Implement iterative/implicit version if needed
        %     [X_adv, Z_adv, U_adv] = attack.dgsm_ev(system_data, true, eps_att);
        
    otherwise
        error('attack:unknownMethod','Unknown attack method: %s', string(method));
end
end
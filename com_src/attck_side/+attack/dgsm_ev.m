function [X_adv, Z_adv, U_adv] = dgsm_ev(sd, use_ev, eps_att)
% dgsm_ev: Direct Gradient Sign Method with Eigenvalue objective
%
% Inputs:
%   sd: SystemData object
%   use_ev: (optional) If true, use eigenvalue-based objective (default: true)
%   eps_att: (optional) Attack parameter (default: from cfg.Const.ATTACKER_UPPERLIMIT)
%
% Outputs:
%   X_adv, Z_adv, U_adv: Adversarial data

if nargin < 2 || isempty(use_ev)
    use_ev = true;
end
if nargin < 3
    eps_att = [];
end

[gradX_pi, gradZ_pi, gradU_pi] = attack.calc_grad(sd);
[X_adv, Z_adv, U_adv] = attack.make_data_adv(sd, gradX_pi, gradZ_pi, gradU_pi);
end

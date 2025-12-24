function [X_adv, Z_adv, U_adv, history] = idgsm_rho_small(sd, eps_att, max_iter, save_history)
% idgsm_rho_small: Iterative Direct Gradient Sign Method to minimize Delta (rho)
%                  using implicit differentiation
%
% This attack minimizes Delta (and thus rho) by using the negative gradient.
% Delta and rho have a monotonic relationship: smaller Delta -> smaller rho.
% 勾配と逆方向にノイズを加える（勾配の符号を反転）。
%
% Inputs:
%   sd: SystemData object
%   eps_att: (optional) Attack parameter (default: from cfg.Const.ATTACKER_UPPERLIMIT)
%   max_iter: (optional) Maximum number of iterations (default: 10)
%   save_history: (optional) If true, save attack history (default: false)
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

% minimize_delta=true を指定してidgsm_implicitを呼び出す（勾配と逆方向）
[X_adv, Z_adv, U_adv, history] = attack.idgsm_implicit(sd, eps_att, max_iter, save_history, true);

end


function [X_adv, Z_adv, U_adv] = make_data_adv(sd, X_grad, Z_grad, U_grad, eps_att, direction)
% make_data_adv: 勾配から攻撃データを生成
%   sd: SystemData object
%   X_grad, Z_grad, U_grad: Delta勾配（各データに対する）
%   eps_att: (optional) 攻撃強度（指定しない場合はcfg.Const.ATTACKER_UPPERLIMITを使用）
%   direction: (optional) 攻撃方向 'positive' (deltaを大きくする) または 'negative' (deltaを小さくする)
%              デフォルトは 'negative'

if nargin < 5 || isempty(eps_att)
    eps_att = cfg.Const.ATTACKER_UPPERLIMIT;
end
if nargin < 6 || isempty(direction)
    direction = 'negative';
end

% 方向の決定
if strcmpi(direction, 'positive')
    direction_mult = 1;  % deltaを大きくする方向（勾配の正方向）
elseif strcmpi(direction, 'negative')
    direction_mult = -1;  % deltaを小さくする方向（勾配の負方向）
else
    error('attack:invalidDirection', 'direction must be ''positive'' or ''negative''');
end

X_adv = d_adv_add(sd.X, X_grad, eps_att, direction_mult);
Z_adv = d_adv_add(sd.Z, Z_grad, eps_att, direction_mult);
U_adv = d_adv_add(sd.U, U_grad, eps_att, direction_mult);

    function d_adv = d_adv_add(d_ori, d_grad, eps, dir_mult)
        % 要素単位で： d_adv = d_ori + ε * sign(d_grad) * dir_mult
        d_adv = d_ori + eps .* sign(d_grad) .* dir_mult;
    end
end
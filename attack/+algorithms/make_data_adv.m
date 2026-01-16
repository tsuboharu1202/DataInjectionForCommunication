function [X_adv, Z_adv, U_adv] = make_data_adv(sd, X_grad, Z_grad, U_grad, eps_att, direction)
% make_data_adv: 勾配から攻撃データを生成
%
% 引数（全て必須）:
%   sd: SystemData オブジェクト
%   X_grad, Z_grad, U_grad: Delta勾配（各データに対する）
%   eps_att: 攻撃強度
%   direction: 攻撃方向 'positive' or 'negative'
%
% 出力:
%   X_adv, Z_adv, U_adv: 攻撃後のデータ

if nargin < 6
    error('make_data_adv:MissingArgs', ...
        '全ての引数（sd, X_grad, Z_grad, U_grad, eps_att, direction）は必須です。');
end

% 方向の決定
if strcmpi(direction, 'positive')
    direction_mult = 1;   % deltaを大きくする方向（勾配の正方向）
elseif strcmpi(direction, 'negative')
    direction_mult = -1;  % deltaを小さくする方向（勾配の負方向）
else
    error('make_data_adv:invalidDirection', ...
        'direction は ''positive'' または ''negative'' を指定してください。');
end

X_adv = d_adv_add(sd.X, X_grad, eps_att, direction_mult);
Z_adv = d_adv_add(sd.Z, Z_grad, eps_att, direction_mult);
U_adv = d_adv_add(sd.U, U_grad, eps_att, direction_mult);

    function d_adv = d_adv_add(d_ori, d_grad, eps, dir_mult)
        % 要素単位で： d_adv = d_ori + ε * sign(d_grad) * dir_mult
        d_adv = d_ori + eps .* sign(d_grad) .* dir_mult;
    end
end

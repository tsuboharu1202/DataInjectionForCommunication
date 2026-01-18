% function [X_adv, Z_adv, U_adv] = dgsm_delta(sd, eps_att, direction, opts)
% % dgsm_delta: Direct Gradient Sign Method with Delta objective
% %
% % Inputs:
% %   sd: SystemData object
% %   eps_att: (optional) 攻撃強度（デフォルト: cfg.Const.ATTACKER_UPPERLIMIT）
% %   direction: (optional) 'positive' (deltaを大きくする) または 'negative' (deltaを小さくする)
% %              デフォルト: 'negative'
% %   opts: (optional) 追加オプション（gammaなど）
% %
% % Outputs:
% %   X_adv, Z_adv, U_adv: Adversarial data

% if nargin < 2
%     eps_att = [];
% end
% if nargin < 3 || isempty(direction)
%     direction = 'negative';
% end
% if nargin < 4
%     opts = struct();
% end

% % Delta勾配を計算（regularizationのみ）
% [gradX_delta, gradZ_delta, gradU_delta] = attack.calc_grad(sd, opts);

% % 攻撃データを生成
% [X_adv, Z_adv, U_adv] = attack.make_data_adv(sd, gradX_delta, gradZ_delta, gradU_delta, eps_att, direction);
% end


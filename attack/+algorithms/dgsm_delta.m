function [X_adv, Z_adv, U_adv] = dgsm_delta(sd, opts)
% dgsm_delta: Direct Gradient Sign Method with Delta objective
%
% 引数:
%   sd: SystemData オブジェクト
%   opts: オプション構造体（全て必須）
%
% 必須オプション:
%   opts.gamma: 正則化パラメータ（SDP計算用）
%   opts.epsilon: 攻撃強度
%   opts.direction: 'positive' or 'negative'
%
% 出力:
%   X_adv, Z_adv, U_adv: 攻撃後のデータ

if nargin < 2
    error('dgsm_delta:MissingOpts', 'opts は必須です。');
end

% 必須パラメータのチェック
required_fields = {'gamma', 'epsilon', 'direction'};
for i = 1:length(required_fields)
    if ~isfield(opts, required_fields{i})
        error('dgsm_delta:MissingField', ...
            'opts.%s は必須です。', required_fields{i});
    end
end

% 勾配を計算
grad_opts = struct('gamma', opts.gamma);
[gradX_delta, gradZ_delta, gradU_delta] = common.calc_grad(sd, grad_opts);

% 攻撃データを生成
[X_adv, Z_adv, U_adv] = algorithms.make_data_adv(sd, gradX_delta, gradZ_delta, gradU_delta, opts.epsilon, opts.direction);
end

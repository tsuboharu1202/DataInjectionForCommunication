function [X_adv, Z_adv, U_adv, history] = execute_attack(system_data, method, opts)
% execute_attack: 攻撃手法を実行
%
% 引数:
%   system_data: SystemData オブジェクト
%   method: 攻撃手法（cfg.AttackType）
%   opts: オプション構造体（全て必須）
%
% 必須オプション（共通）:
%   opts.gamma: 正則化パラメータ（SDP計算用）
%   opts.epsilon: 攻撃強度
%   opts.direction: 'positive' or 'negative'
%
% IDGSM用の追加必須オプション:
%   opts.alpha: ステップサイズ
%   opts.max_iteration: 最大反復回数
%
% オプション:
%   opts.save_history: 履歴を保存するか（デフォルト: false）
%   opts.save_grad_norms: 勾配ノルムを保存するか（デフォルト: false）
%   opts.normalize_grad: 勾配を正規化するか（デフォルト: true）
%
% 出力:
%   X_adv, Z_adv, U_adv: 攻撃後のデータ
%   history: 履歴（IDGSM + save_history=true の場合のみ）

if nargin < 3
    error('execute_attack:MissingOpts', 'opts は必須です。');
end

% 共通の必須パラメータのチェック
required_fields = {'gamma', 'epsilon', 'direction'};
for i = 1:length(required_fields)
    if ~isfield(opts, required_fields{i})
        error('execute_attack:MissingField', ...
            'opts.%s は必須です。', required_fields{i});
    end
end

history = [];

% methodを文字列に変換
method_str = string(method);

switch method_str
    case {cfg.AttackType.DIRECT_DGSM_DELTA, "DIRECT_DGSM_DELTA"}
        [X_adv, Z_adv, U_adv] = algorithms.dgsm_delta(system_data, opts);
        
    case {cfg.AttackType.IMPLICIT_IDGSM_DELTA, "IMPLICIT_IDGSM_DELTA"}
        % IDGSM用の追加必須パラメータをチェック
        idgsm_required = {'alpha', 'max_iteration'};
        for i = 1:length(idgsm_required)
            if ~isfield(opts, idgsm_required{i})
                error('execute_attack:MissingField', ...
                    'IDGSM では opts.%s は必須です。', idgsm_required{i});
            end
        end
        [X_adv, Z_adv, U_adv, history] = algorithms.idgsm_delta(system_data, opts);
        
    otherwise
        error('attack:unknownMethod', '不明な攻撃手法です: %s', method_str);
end
end

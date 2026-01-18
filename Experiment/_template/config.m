% config.m - 実験設定ファイル
% このファイルを編集して実験パラメータを設定してください
%
% 使い方:
%   1. このテンプレートフォルダをコピーして新しい実験フォルダを作成
%   2. このconfig.mを編集
%   3. run_experiment.mを実行

%% ========== 実験情報 ==========
cfg = struct();
cfg.experiment_name = 'Experiment_Template';  % 実験名
cfg.description = '実験の説明をここに書く';   % 説明
cfg.date = datestr(now, 'yyyy-mm-dd');        % 実験日

%% ========== システム設定 ==========
% システム行列（編集してください）
cfg.A = [-0.192, -0.936, -0.814;
         -0.918,  0.729, -0.724;
         -0.412, -0.135, -0.516];
cfg.B = [-0.554; 0.735; 0.528];

% システムサイズ（自動計算）
cfg.n = size(cfg.A, 1);
cfg.m = size(cfg.B, 2);

%% ========== データ生成設定 ==========
cfg.T = 20;  % サンプル数（cfg.Const.SAMPLE_COUNTを使う場合はコメントアウト）

%% ========== 実験パラメータ ==========
% 比較するパラメータ（実験に応じて変更）
cfg.param_list = [1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1];

% 攻撃設定
cfg.attack_eps = 5e-3;      % 攻撃の大きさ
cfg.attack_direction = 'positive';  % 'positive' or 'negative'
cfg.gamma_attack = 1e3;     % 攻撃計算時のgamma

% SDPパラメータ
cfg.gamma_sdp = 1e3;        % SDPのgammaパラメータ

% 試行回数
cfg.num_trials = 10;

%% ========== 手法設定 ==========
% 比較する手法（'baseline', 'proposed', 'proposed_robust'）
cfg.methods = {'baseline', 'proposed'};

%% ========== 保存設定 ==========
cfg.save_results = true;    % 結果を保存するか
cfg.save_figures = true;    % 図を保存するか
cfg.figure_format = 'png';  % 図のフォーマット ('png', 'fig', 'pdf')

%% ========== 判定条件 ==========
cfg.tolerance = 1.01;       % 判定時のマージン（1.01 = 1%マージン）

%% ========== その他 ==========
cfg.verbose = true;         % 詳細出力
cfg.random_seed = [];       % 乱数シード（[]で毎回異なる、数値で固定）


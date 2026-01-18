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
cfg.description = '同一のシステム(cfg)を用いて、次回以降のシミュレーションに使うための攻撃データの生成';   % 説明

cfg.date = datestr(now, '2026-01-18');        % 実験日

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
cfg.T = 10;  % サンプル数（cfg.Const.SAMPLE_COUNTを使う場合はコメントアウト）

%% ========== SDP設定 ==========
cfg.gamma = 1e3;  % 正則化パラメータ（固定）

%% ========== 実験パラメータ ==========
% 比較するパラメータ（実験に応じて変更）
cfg.attack_eps_list = [1e-4, 1e-3, 1e-2, 1e-1, 1e0];
cfg.step_size_coefficient = 1e-1;
cfg.step_max_iter = 30;
cfg.trial = 20;

%% ========== 手法設定 ==========
% 比較する手法（'baseline', 'proposed', 'proposed_robust'）
% cfg.methods = {'baseline', 'proposed'};

%% ========== 保存設定 ==========
cfg.save_results = false;    % 結果を保存するか
cfg.save_figures = false;    % 図を保存するか
cfg.figure_format = 'png';  % 図のフォーマット ('png', 'fig', 'pdf')

%% ========== 判定条件 ==========
% cfg.tolerance = 1.01;       % 判定時のマージン（1.01 = 1%マージン）

%% ========== その他 ==========
cfg.verbose = true;         % 詳細出力
cfg.random_seed = [];       % 乱数シード（[]で毎回異なる、数値で固定）


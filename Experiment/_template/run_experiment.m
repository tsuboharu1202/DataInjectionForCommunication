% run_experiment.m - 実験実行テンプレート
%
% 使い方:
%   1. config.m を編集
%   2. このスクリプトを実行
%   3. 結果は results/ フォルダに保存される
%
clear; clc; close all;

%% ========== 設定読み込み ==========
run('config.m');

fprintf('==============================================\n');
fprintf('  実験: %s\n', cfg.experiment_name);
fprintf('  日時: %s\n', datestr(now));
fprintf('==============================================\n\n');

%% ========== 乱数シード設定 ==========
if ~isempty(cfg.random_seed)
    rng(cfg.random_seed);
    fprintf('乱数シード: %d\n', cfg.random_seed);
end

%% ========== 結果フォルダ確認 ==========
if ~exist('results', 'dir'), mkdir('results'); end
if ~exist('results/figures', 'dir'), mkdir('results/figures'); end

%% ========== ログ開始 ==========
log_filename = sprintf('results/log_%s.txt', datestr(now, 'yyyymmdd_HHMMSS'));
diary(log_filename);
fprintf('ログファイル: %s\n\n', log_filename);

%% ========== 設定の表示 ==========
if cfg.verbose
    fprintf('--- 実験設定 ---\n');
    fprintf('システムサイズ: n=%d, m=%d, T=%d\n', cfg.n, cfg.m, cfg.T);
    fprintf('パラメータ: ');
    fprintf('%.0e ', cfg.param_list);
    fprintf('\n');
    fprintf('攻撃サイズ: %.0e\n', cfg.attack_eps);
    fprintf('試行回数: %d\n', cfg.num_trials);
    fprintf('手法: ');
    fprintf('%s ', cfg.methods{:});
    fprintf('\n\n');
end

%% ========== rho_rough計算 ==========
eigA = eig(cfg.A);
gamma_prod = 1;
for i = 1:numel(eigA)
    if abs(eigA(i)) > 1
        gamma_prod = gamma_prod * abs(eigA(i));
    end
end
rho_rough = 1 / gamma_prod;
fprintf('rho_rough (理論的下界): %.6e\n\n', rho_rough);

%% ========== 結果格納 ==========
results = struct();
results.cfg = cfg;
results.rho_rough = rho_rough;
results.timestamp = datestr(now);

% ここに実験結果を格納する配列を初期化
% 例:
% results.delta = zeros(length(cfg.param_list), cfg.num_trials);
% results.success_rate = zeros(length(cfg.param_list), 1);

%% ========== メイン実験ループ ==========
fprintf('========== 実験開始 ==========\n\n');

% TODO: ここに実験ロジックを実装
% 
% 例:
% for trial = 1:cfg.num_trials
%     fprintf('--- Trial %d/%d ---\n', trial, cfg.num_trials);
%     
%     % データ生成
%     V = input.make_inputU(cfg.m);
%     [X, Z] = datasim.simulate_openloop(cfg.A, cfg.B, V);
%     U = V;
%     
%     for p_idx = 1:length(cfg.param_list)
%         param = cfg.param_list(p_idx);
%         
%         % SDP解く
%         % ...
%     end
% end

warning('実験ロジックが未実装です。run_experiment.mを編集してください。');

%% ========== 結果サマリ ==========
fprintf('\n========== 結果サマリ ==========\n');
% TODO: 結果のサマリを表示

%% ========== 結果保存 ==========
if cfg.save_results
    results_filename = 'results/data.mat';
    save(results_filename, 'results', 'cfg', 'rho_rough');
    fprintf('\n結果を保存しました: %s\n', results_filename);
end

%% ========== 描画 ==========
if cfg.save_figures
    fprintf('\n描画中...\n');
    run('plot_results.m');
end

%% ========== 終了 ==========
diary off;
fprintf('\n==============================================\n');
fprintf('  実験完了: %s\n', datestr(now));
fprintf('==============================================\n');


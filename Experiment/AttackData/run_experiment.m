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
    fprintf('正則化パラメータ: gamma=%.0e\n', cfg.gamma);
    fprintf('攻撃サイズ (epsilon): ');
    fprintf('%.0e ', cfg.attack_eps_list);
    fprintf('\n');
    fprintf('ステップサイズ係数: %.0e\n', cfg.step_size_coefficient);
    fprintf('最大反復回数: %d\n', cfg.step_max_iter);
    fprintf('試行回数: %d\n', cfg.trial);
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

% Phi行列の設定（IQC用）
Phi11 = 1e-5 * eye(cfg.n);
Phi12 = zeros(cfg.n, cfg.T);
Phi22 = -eye(cfg.T);

% 結果格納用の構造体を初期化
results.datasets = struct([]);

% 20セットのデータ生成と攻撃実行
for trial = 1:cfg.trial
    fprintf('========================================\n');
    fprintf('--- Trial %d/%d ---\n', trial, cfg.trial);
    fprintf('========================================\n');
    
    % ========== データ生成 ==========
    V = core.make_inputU(cfg.m, cfg.T);
    [X, Z, U] = datasim.simulate_openloop_stable(cfg.A, cfg.B, V);
    
    % SystemData作成
    data_original = datasim.SystemData(cfg.A, cfg.B, X, Z, U, Phi11, Phi12, Phi22);
    
    % 元のデータでSDPを解く
    try
        [sol_original, ~, ~, ~, diag_original] = proposed.solve_sdp(data_original, cfg.gamma);
        if diag_original.problem ~= 0
            warning('Trial %d: 元のデータでSDPが解けませんでした (status=%d)', trial, diag_original.problem);
            continue;
        end
        delta_original = sol_original.delta;
    catch ME
        warning('Trial %d: SDPエラー - %s', trial, ME.message);
        continue;
    end
    
    fprintf('  元のdelta: %.6e\n', delta_original);
    
    % 結果構造体に元のデータを保存
    results.datasets(trial).original = struct();
    results.datasets(trial).original.U = U;
    results.datasets(trial).original.X = X;
    results.datasets(trial).original.Z = Z;
    results.datasets(trial).original.delta = delta_original;
    results.datasets(trial).original.K = sol_original.K;
    
    % 攻撃結果を格納する配列
    results.datasets(trial).attacks = struct([]);
    attack_idx = 0;
    
    % ========== 各epsilonで攻撃実行（negative と positive 両方） ==========
    for eps_idx = 1:length(cfg.attack_eps_list)
        epsilon = cfg.attack_eps_list(eps_idx);
        step_size = epsilon * cfg.step_size_coefficient;
        
        fprintf('\n  --- Epsilon = %.0e (step_size=%.0e, max_iter=%d) ---\n', ...
            epsilon, step_size, cfg.step_max_iter);
        
        % negative と positive の両方を実行
        for dir_idx = 1:2
            if dir_idx == 1
                direction = 'negative';
            else
                direction = 'positive';
            end
            
            fprintf('    Direction: %s\n', direction);
            
            % IDGSM攻撃オプション設定
            opts_attack = struct();
            opts_attack.gamma = cfg.gamma;
            opts_attack.epsilon = epsilon;
            opts_attack.alpha = step_size;
            opts_attack.max_iteration = cfg.step_max_iter;
            opts_attack.direction = direction;
            opts_attack.save_history = false;  % 履歴は保存しない（データが大きいため）
            opts_attack.save_grad_norms = false;
            opts_attack.normalize_grad = false;  % 正規化しない
            
            % 攻撃実行
            try
                [X_adv, Z_adv, U_adv, history] = algorithms.idgsm_delta(data_original, opts_attack);
                
                % 攻撃後のデータでSDPを解く
                data_attacked = datasim.SystemData(cfg.A, cfg.B, X_adv, Z_adv, U_adv, Phi11, Phi12, Phi22);
                [sol_attacked, ~, ~, ~, diag_attacked] = proposed.solve_sdp(data_attacked, cfg.gamma);
                
                if diag_attacked.problem ~= 0
                    warning('      Epsilon=%.0e, %s: 攻撃後のSDPが解けませんでした (status=%d)', ...
                        epsilon, direction, diag_attacked.problem);
                    delta_attacked = NaN;
                    K_attacked = NaN;
                else
                    delta_attacked = sol_attacked.delta;
                    K_attacked = sol_attacked.K;
                end
                
                fprintf('      攻撃後delta: %.6e (変化: %.2f%%)\n', ...
                    delta_attacked, 100*(delta_attacked - delta_original)/delta_original);
                
                % 結果を保存
                attack_idx = attack_idx + 1;
                results.datasets(trial).attacks(attack_idx).epsilon = epsilon;
                results.datasets(trial).attacks(attack_idx).step_size = step_size;
                results.datasets(trial).attacks(attack_idx).direction = direction;
                results.datasets(trial).attacks(attack_idx).U_adv = U_adv;
                results.datasets(trial).attacks(attack_idx).X_adv = X_adv;
                results.datasets(trial).attacks(attack_idx).Z_adv = Z_adv;
                results.datasets(trial).attacks(attack_idx).delta = delta_attacked;
                results.datasets(trial).attacks(attack_idx).K = K_attacked;
                results.datasets(trial).attacks(attack_idx).success = (diag_attacked.problem == 0);
                
            catch ME
                warning('      Epsilon=%.0e, %s: 攻撃実行エラー - %s', epsilon, direction, ME.message);
                % エラー時はデータを記録しない（スキップ）
                continue;
                
            end
        end
    end
    
    % ========== 各データセットごとに保存 ==========
    if cfg.save_results
        dataset_filename = sprintf('results/dataset_%03d_%s.mat', trial, datestr(now, 'yyyymmdd_HHMMSS'));
        dataset_data = struct();
        dataset_data.original = results.datasets(trial).original;
        dataset_data.attacks = results.datasets(trial).attacks;
        dataset_data.cfg = cfg;
        dataset_data.rho_rough = rho_rough;
        dataset_data.trial = trial;
        dataset_data.timestamp = datestr(now);
        
        save(dataset_filename, 'dataset_data', '-v7.3');
        fprintf('  データセット %d を保存しました: %s\n', trial, dataset_filename);
    end
    
    fprintf('\n');
end

%% ========== 結果サマリ ==========
fprintf('\n========== 結果サマリ ==========\n');

% 成功した攻撃の統計
success_count = 0;
total_count = 0;
delta_changes = [];

for trial = 1:cfg.trial
    if isfield(results.datasets, 'attacks') && ~isempty(results.datasets(trial).attacks)
        for eps_idx = 1:length(results.datasets(trial).attacks)
            if results.datasets(trial).attacks(eps_idx).success
                success_count = success_count + 1;
                delta_orig = results.datasets(trial).original.delta;
                delta_adv = results.datasets(trial).attacks(eps_idx).delta;
                if ~isnan(delta_adv) && ~isnan(delta_orig)
                    delta_changes(end+1) = (delta_adv - delta_orig) / delta_orig * 100;
                end
            end
            total_count = total_count + 1;
        end
    end
end

fprintf('総攻撃試行数: %d\n', total_count);
fprintf('成功数: %d (成功率: %.1f%%)\n', success_count, 100*success_count/total_count);
if ~isempty(delta_changes)
    fprintf('Delta変化の統計:\n');
    fprintf('  平均: %.2f%%\n', mean(delta_changes));
    fprintf('  最小: %.2f%%\n', min(delta_changes));
    fprintf('  最大: %.2f%%\n', max(delta_changes));
end

%% ========== 結果保存（まとめて保存する場合） ==========
% 各データセットは個別に保存済み
% 必要に応じて、全データをまとめて保存する場合は以下を有効化
% if cfg.save_results
%     results_filename = sprintf('results/attack_data_all_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
%     save(results_filename, 'results', 'cfg', 'rho_rough', '-v7.3');
%     fprintf('\n全データをまとめて保存しました: %s\n', results_filename);
% end

%% ========== 終了 ==========
diary off;
fprintf('\n==============================================\n');
fprintf('  実験完了: %s\n', datestr(now));
fprintf('==============================================\n');


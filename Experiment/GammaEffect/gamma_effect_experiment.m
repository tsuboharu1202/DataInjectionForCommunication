% gamma_effect_experiment.m
% Gammaパラメータが攻撃成功率に与える影響を評価する実験
%
% 実験内容:
% - システム: demo_implicit_regularization.m (19-27) のA, B
% - 攻撃の大きさ: 1e-5, 1e-4, 1e-3（それぞれ10回ずつ）
% - IDGSMで30 iter、positive方向（deltaを大きくする）
% - stepsize = epsilon / 5
% - gamma = 0, 1e3, 1e5 の3通りで評価
% - 攻撃成功判定: 1/delta < hinfnorm_AK(A, B, K)

clear; clc; close all;

fprintf('=== Gamma効果評価実験 ===\n\n');

% ============================================
% 1. システム定義
% ============================================
fprintf('1. システム定義...\n');
A = [-0.192, -0.936, -0.814;
    -0.918,  0.729, -0.724;
    -0.412,  0.735, -0.516];
B = [-0.554;
    0.735;
    0.528];

n = size(A, 1);  % 状態次元
m = size(B, 2);  % 入力次元（Bの列数）
T = 20;  % サンプル数

fprintf('  システムサイズ: n=%d, m=%d, T=%d\n', n, m, T);
fprintf('  Aの固有値:\n');
disp(eig(A));

% ============================================
% 2. 実験パラメータ
% ============================================
fprintf('\n2. 実験パラメータ設定...\n');

% 攻撃の大きさ（epsilon）
epsilon_list = [1e-5, 1e-4, 1e-3];
n_epsilon = length(epsilon_list);
n_trials_per_epsilon = 10;  % 各epsilonについて10回ずつ

% Gamma値
gamma_list = [1e-4,1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5];
n_gamma = length(gamma_list);

% IDGSMパラメータ
max_iter = 30;  % 最大反復回数
direction = 'positive';  % deltaを大きくする方向
max_data_attempts = 20;  % データ生成の最大試行回数（infeasible対策）

% Phi設定
Phi11 = 1e-1 * eye(n);
Phi12 = zeros(n, T);
Phi22 = -eye(T);

fprintf('  攻撃の大きさ: ');
fprintf('%.0e ', epsilon_list);
fprintf('\n');
fprintf('  各epsilonの試行回数: %d\n', n_trials_per_epsilon);
fprintf('  Gamma値: ');
fprintf('%.0e ', gamma_list);
fprintf('\n');
fprintf('  IDGSM最大反復回数: %d\n', max_iter);
fprintf('  攻撃方向: %s\n', direction);

% ============================================
% 3. 結果保存用の構造体配列
% ============================================
results = struct();
results.epsilon = [];
results.trial_idx = [];
results.gamma = [];
results.delta = [];
results.hinf_norm = [];
results.attack_success = [];  % 1/delta < hinf_norm
results.infeasible = [];  % SDPがinfeasibleになったかどうか

% ============================================
% 4. 実験実行
% ============================================
fprintf('\n3. 実験実行開始...\n');
total_trials = n_epsilon * n_trials_per_epsilon;
trial_count = 0;

for eps_idx = 1:n_epsilon
    epsilon = epsilon_list(eps_idx);
    stepsize =  epsilon*10; % ステップサイズ = epsilon / 5
    
    fprintf('\n--- Epsilon = %.0e ---\n', epsilon);
    
    for trial_idx = 1:n_trials_per_epsilon
        trial_count = trial_count + 1;
        fprintf('\n[試行 %d/%d] Epsilon=%.0e, Trial=%d/%d\n', ...
            trial_count, total_trials, epsilon, trial_idx, n_trials_per_epsilon);
        
        % ============================================
        % 4.1 データ生成（infeasible対策）
        % ============================================
        data_valid = false;
        data_attempt = 0;
        data_ori = [];
        
        while ~data_valid && data_attempt < max_data_attempts
            data_attempt = data_attempt + 1;
            
            try
                % データ生成
                V = make_inputU(m, T);
                [X, Z, U] = datasim.simulate_openloop_stable(A, B, V);
                
                data_ori = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
                
                % 初期SDPを解いてデータが有効か確認（gamma=0で）
                [sol_init, ~, ~, ~, diagnostics_init] = proposed.solve_sdp(data_ori, 1);
                
                if diagnostics_init.problem == 0
                    data_valid = true;
                    fprintf('  データ生成成功 (試行 %d)\n', data_attempt);
                else
                    fprintf('  データ生成試行 %d: SDP infeasible (status=%d)\n', ...
                        data_attempt, diagnostics_init.problem);
                end
            catch ME
                fprintf('  データ生成試行 %d: エラー - %s\n', data_attempt, ME.message);
            end
        end
        
        if ~data_valid
            fprintf('  警告: データ生成に失敗しました（最大試行回数に達しました）\n');
            % この試行をスキップ
            for g_idx = 1:n_gamma
                gamma = gamma_list(g_idx);
                results.epsilon(end+1) = epsilon;
                results.trial_idx(end+1) = trial_idx;
                results.gamma(end+1) = gamma;
                results.delta(end+1) = NaN;
                results.hinf_norm(end+1) = NaN;
                results.attack_success(end+1) = false;
                results.infeasible(end+1) = true;
            end
            continue;
        end
        
        % ============================================
        % 4.2 IDGSM攻撃を実行（gamma=1e3で数値的に安定に実行）
        % ============================================
        fprintf('  IDGSM攻撃を実行（gamma=1）...\n');
        
        % 攻撃オプション
        opts_attack = struct();
        opts_attack.gamma = 1;
        opts_attack.epsilon = epsilon;
        opts_attack.alpha = stepsize;
        opts_attack.direction = direction;
        opts_attack.save_history = true;
        
        try
            [X_attacked, Z_attacked, U_attacked, history_attack] = algorithms.idgsm_delta(data_ori, opts_attack);
            
            % history_attackが空でない場合のみ反復回数を表示
            if ~isempty(history_attack) && isstruct(history_attack) && isfield(history_attack, 'iter_count')
                fprintf('  攻撃完了 (反復回数: %d)\n', history_attack.iter_count);
            else
                fprintf('  攻撃完了\n');
            end
            
            % 攻撃後のデータ
            data_attacked = datasim.SystemData(A, B, X_attacked, Z_attacked, U_attacked, ...
                Phi11, Phi12, Phi22);
            
        catch ME
            fprintf('  警告: IDGSM攻撃中にエラー - %s\n', ME.message);
            % この試行をスキップ
            for g_idx = 1:n_gamma
                gamma = gamma_list(g_idx);
                results.epsilon(end+1) = epsilon;
                results.trial_idx(end+1) = trial_idx;
                results.gamma(end+1) = gamma;
                results.delta(end+1) = NaN;
                results.hinf_norm(end+1) = NaN;
                results.attack_success(end+1) = false;
                results.infeasible(end+1) = true;
            end
            continue;
        end
        
        % ============================================
        % 4.3 攻撃後のデータでgammaを3通り変えてSDPを解く
        % ============================================
        fprintf('  攻撃後のデータでgammaを変えてSDPを解く...\n');
        
        for g_idx = 1:n_gamma
            gamma = gamma_list(g_idx);
            fprintf('    gamma = %.0e: ', gamma);
            
            try
                [sol, K, ~, ~, diagnostics] = proposed.solve_sdp(data_attacked, gamma);
                
                if diagnostics.problem == 0
                    delta = sol.delta;
                    hinf_norm = helper.hinfnorm_AK(A, B, K);
                    attack_success = (1/delta < hinf_norm);
                    infeasible = false;
                    
                    fprintf('delta=%.6e, hinf_norm=%.6e, success=%d\n', ...
                        delta, hinf_norm, attack_success);
                else
                    fprintf('SDP infeasible (status=%d)\n', diagnostics.problem);
                    delta = NaN;
                    hinf_norm = NaN;
                    attack_success = false;
                    infeasible = true;
                end
            catch ME
                fprintf('エラー - %s\n', ME.message);
                delta = NaN;
                hinf_norm = NaN;
                attack_success = false;
                infeasible = true;
            end
            
            % 結果を保存
            results.epsilon(end+1) = epsilon;
            results.trial_idx(end+1) = trial_idx;
            results.gamma(end+1) = gamma;
            results.delta(end+1) = delta;
            results.hinf_norm(end+1) = hinf_norm;
            results.attack_success(end+1) = attack_success;
            results.infeasible(end+1) = infeasible;
        end
    end
end

% ============================================
% 5. 結果の集計と表示
% ============================================
fprintf('\n\n=== 実験結果 ===\n\n');

% 結果をテーブルに変換（解析しやすくするため）
results_table = struct2table(results);

% 各epsilon、各gammaについて成功率を計算
fprintf('攻撃成功率（1/delta < hinf_norm）:\n');
fprintf('%-10s | %-10s | %-15s | %-15s | %-15s\n', ...
    'Epsilon', 'Gamma', '成功率', '有効試行数', 'Infeasible数');
fprintf('%s\n', repmat('-', 1, 80));

for eps_idx = 1:n_epsilon
    epsilon = epsilon_list(eps_idx);
    
    for g_idx = 1:n_gamma
        gamma = gamma_list(g_idx);
        
        % 該当する試行を抽出
        mask = (results_table.epsilon == epsilon) & (results_table.gamma == gamma);
        subset = results_table(mask, :);
        
        % 有効な試行（infeasibleでない）を抽出
        valid_mask = ~subset.infeasible;
        valid_subset = subset(valid_mask, :);
        
        n_valid = sum(valid_mask);
        n_infeasible = sum(subset.infeasible);
        n_success = sum(valid_subset.attack_success);
        
        if n_valid > 0
            success_rate = n_success / n_valid * 100;
        else
            success_rate = NaN;
        end
        
        fprintf('%-10.0e | %-10.0e | %-13.1f%% | %-15d | %-15d\n', ...
            epsilon, gamma, success_rate, n_valid, n_infeasible);
    end
end

% ============================================
% 6. 結果を保存
% ============================================
fprintf('\n結果を保存中...\n');
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
save_filename = fullfile('Experiment', sprintf('gamma_effect_results_%s.mat', timestamp));
save(save_filename, 'results', 'results_table', 'epsilon_list', 'gamma_list', ...
    'A', 'B', 'n', 'm', 'T', 'max_iter', 'direction');
fprintf('  保存先: %s\n', save_filename);

% CSV形式でも保存（解析しやすくするため）
csv_filename = fullfile('Experiment', sprintf('gamma_effect_results_%s.csv', timestamp));
writetable(results_table, csv_filename);
fprintf('  CSV保存先: %s\n', csv_filename);

fprintf('\n実験完了！\n');


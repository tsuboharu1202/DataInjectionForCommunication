% run_experiment.m
% original_thesis と regularization_sdp の性能比較実験
%
% 実験1: 保守性の評価（攻撃なし）
% 実験2: ロバスト性の評価（攻撃あり）
%
clear; clc; close all;

%% ========== パラメータ設定 ==========
% システム設定（cfg.Systemから取得）
A = cfg.System.A;
B = cfg.System.B;
[n, m] = cfg.System.getDimensions();
T = cfg.System.getSampleCount();

% パラメータ範囲
param_list = [1e-4, 1e-2, 1];  % Phi11係数 / epsilon

% 実験1の設定
num_trials_exp1 = 10;  % 10回の平均

% 実験2の設定
attack_sizes = [1e-3, 1e-2, 1e-1];  % 攻撃の大きさ
num_trials_exp2 = 5;  % 5回の平均
gamma_attack = 1e3;  % IDGSM攻撃時のgamma
gamma_reg = 1e3;     % regularization_sdpのgamma（固定）

% rho_roughの計算（cfg.Systemのヘルパー関数を使用）
rho_rough = cfg.System.calcRhoRough(A);
fprintf('rho_rough (理論的下界): %.6e\n', rho_rough);

%% ========== 実験1: 保守性の評価（攻撃なし） ==========
fprintf('\n========== 実験1: 保守性の評価 ==========\n');

% 結果格納
results_exp1 = struct();
results_exp1.param_list = param_list;
results_exp1.delta_ori = zeros(length(param_list), num_trials_exp1);  % original_thesis
results_exp1.delta_reg = zeros(length(param_list), num_trials_exp1);  % regularization_sdp
results_exp1.infeasible_ori = false(length(param_list), num_trials_exp1);
results_exp1.infeasible_reg = false(length(param_list), num_trials_exp1);

for trial = 1:num_trials_exp1
    fprintf('\n--- Trial %d/%d ---\n', trial, num_trials_exp1);
    
    % データ生成
    V = make_inputU(m, T);
    [X, Z] = datasim.simulate_openloop(A, B, V);
    U = V;
    
    for p_idx = 1:length(param_list)
        param = param_list(p_idx);
        fprintf('  param = %.0e: ', param);
        
        % Phi行列の設定（original_thesis用）
        Phi11 = param * eye(n);
        Phi12 = zeros(n, T);
        Phi22 = -eye(T);
        
        % データ構造作成
        data = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
        
        % --- original_thesis ---
        try
            [sol_ori, K_ori, ~, ~, diag_ori] = baseline.solve_sdp(data);
            results_exp1.delta_ori(p_idx, trial) = sol_ori.delta;
            fprintf('ori=%.4e, ', sol_ori.delta);
        catch ME
            fprintf('ori=INFEASIBLE, ');
            results_exp1.delta_ori(p_idx, trial) = NaN;
            results_exp1.infeasible_ori(p_idx, trial) = true;
        end
        
        % --- regularization_sdp (with robust) ---
        try
            [sol_reg, K_reg, delta_reg, ~, diag_reg] = proposed.solve_sdp_with_robust(data, param, gamma_reg);
            results_exp1.delta_reg(p_idx, trial) = delta_reg;
            fprintf('reg=%.4e\n', delta_reg);
        catch ME
            fprintf('reg=INFEASIBLE\n');
            results_exp1.delta_reg(p_idx, trial) = NaN;
            results_exp1.infeasible_reg(p_idx, trial) = true;
        end
    end
end

% 実験1の平均を計算
results_exp1.delta_ori_mean = nanmean(results_exp1.delta_ori, 2);
results_exp1.delta_reg_mean = nanmean(results_exp1.delta_reg, 2);
results_exp1.delta_ori_std = nanstd(results_exp1.delta_ori, 0, 2);
results_exp1.delta_reg_std = nanstd(results_exp1.delta_reg, 0, 2);
results_exp1.rho_rough = rho_rough;

fprintf('\n=== 実験1 結果サマリ ===\n');
fprintf('param\t\toriginal_thesis\t\tregularization_sdp\n');
for p_idx = 1:length(param_list)
    fprintf('%.0e\t\t%.4e (±%.2e)\t%.4e (±%.2e)\n', ...
        param_list(p_idx), ...
        results_exp1.delta_ori_mean(p_idx), results_exp1.delta_ori_std(p_idx), ...
        results_exp1.delta_reg_mean(p_idx), results_exp1.delta_reg_std(p_idx));
end

%% ========== 実験2: ロバスト性の評価（攻撃あり） ==========
fprintf('\n========== 実験2: ロバスト性の評価 ==========\n');

% 結果格納
results_exp2 = struct();
results_exp2.param_list = param_list;
results_exp2.attack_sizes = attack_sizes;
results_exp2.rho_rough = rho_rough;

% 判定結果の格納 [param_idx, attack_idx, trial]
% original_thesis
results_exp2.ori_delta = zeros(length(param_list), length(attack_sizes), num_trials_exp2);
results_exp2.ori_hinf = zeros(length(param_list), length(attack_sizes), num_trials_exp2);
results_exp2.ori_fail_delta = false(length(param_list), length(attack_sizes), num_trials_exp2);
results_exp2.ori_fail_hinf = false(length(param_list), length(attack_sizes), num_trials_exp2);
results_exp2.ori_infeasible = false(length(param_list), length(attack_sizes), num_trials_exp2);

% regularization_sdp
results_exp2.reg_delta = zeros(length(param_list), length(attack_sizes), num_trials_exp2);
results_exp2.reg_hinf = zeros(length(param_list), length(attack_sizes), num_trials_exp2);
results_exp2.reg_fail_delta = false(length(param_list), length(attack_sizes), num_trials_exp2);
results_exp2.reg_fail_hinf = false(length(param_list), length(attack_sizes), num_trials_exp2);
results_exp2.reg_infeasible = false(length(param_list), length(attack_sizes), num_trials_exp2);

for trial = 1:num_trials_exp2
    fprintf('\n===== Trial %d/%d =====\n', trial, num_trials_exp2);
    
    % データ生成（攻撃前）
    V = make_inputU(m, T);
    [X_ori, Z_ori] = datasim.simulate_openloop(A, B, V);
    U_ori = V;
    
    % 基本のPhi行列（攻撃計算用）
    Phi11_base = 1e-1 * eye(n);
    Phi12 = zeros(n, T);
    Phi22 = -eye(T);
    
    for atk_idx = 1:length(attack_sizes)
        attack_eps = attack_sizes(atk_idx);
        fprintf('\n--- Attack size: %.0e ---\n', attack_eps);
        
        % IDGSM攻撃を実行
        data_base = datasim.SystemData(A, B, X_ori, Z_ori, U_ori, Phi11_base, Phi12, Phi22);
        
        % 攻撃オプション
        opts_attack = struct();
        opts_attack.gamma = gamma_attack;
        opts_attack.epsilon = attack_eps;
        opts_attack.alpha = attack_eps / 5;
        opts_attack.direction = 'positive';
        
        try
            [X_atk, Z_atk, U_atk, ~] = algorithms.idgsm_delta(data_base, opts_attack);
            fprintf('  攻撃完了\n');
        catch ME
            % 攻撃失敗時は明示的に表示し、この試行の結果はすべてNaN
            fprintf('  [攻撃失敗] %s\n', ME.message);
            fprintf('  → この攻撃サイズの試行をスキップします（結果はNaN）\n');
            
            % すべてのパラメータに対してNaNを記録
            for p_idx = 1:length(param_list)
                results_exp2.ori_delta(p_idx, atk_idx, trial) = NaN;
                results_exp2.ori_hinf(p_idx, atk_idx, trial) = NaN;
                results_exp2.ori_infeasible(p_idx, atk_idx, trial) = true;
                results_exp2.reg_delta(p_idx, atk_idx, trial) = NaN;
                results_exp2.reg_hinf(p_idx, atk_idx, trial) = NaN;
                results_exp2.reg_infeasible(p_idx, atk_idx, trial) = true;
            end
            continue;  % 次の攻撃サイズへ
        end
        
        % 各パラメータでSDP評価
        for p_idx = 1:length(param_list)
            param = param_list(p_idx);
            fprintf('  param = %.0e: ', param);
            
            % Phi行列の設定
            Phi11 = param * eye(n);
            data_attacked = datasim.SystemData(A, B, X_atk, Z_atk, U_atk, Phi11, Phi12, Phi22);
            
            % --- original_thesis ---
            try
                [sol_ori, K_ori, ~, ~, ~] = baseline.solve_sdp(data_attacked);
                delta_ori = sol_ori.delta;
                hinf_ori = helper.hinfnorm_AK(A, B, K_ori);
                
                results_exp2.ori_delta(p_idx, atk_idx, trial) = delta_ori;
                results_exp2.ori_hinf(p_idx, atk_idx, trial) = hinf_ori;
                
                % 判定1: delta > rho_rough * 1.01
                if delta_ori > rho_rough * 1.01
                    results_exp2.ori_fail_delta(p_idx, atk_idx, trial) = true;
                end
                % 判定2: hinf_norm > (1/delta) * 1.01
                if hinf_ori > (1/delta_ori) * 1.01
                    results_exp2.ori_fail_hinf(p_idx, atk_idx, trial) = true;
                end
                
                fprintf('ori: delta=%.3e, hinf=%.3e, ', delta_ori, hinf_ori);
            catch ME
                fprintf('ori: INFEASIBLE, ');
                results_exp2.ori_delta(p_idx, atk_idx, trial) = NaN;
                results_exp2.ori_hinf(p_idx, atk_idx, trial) = NaN;
                results_exp2.ori_infeasible(p_idx, atk_idx, trial) = true;
                results_exp2.ori_fail_delta(p_idx, atk_idx, trial) = true;
                results_exp2.ori_fail_hinf(p_idx, atk_idx, trial) = true;
            end
            
            % --- regularization_sdp (with robust) ---
            try
                [sol_reg, K_reg, delta_reg, ~, ~] = proposed.solve_sdp_with_robust(data_attacked, param, gamma_reg);
                hinf_reg = helper.hinfnorm_AK(A, B, K_reg);
                
                results_exp2.reg_delta(p_idx, atk_idx, trial) = delta_reg;
                results_exp2.reg_hinf(p_idx, atk_idx, trial) = hinf_reg;
                
                % 判定1: delta > rho_rough * 1.01
                if delta_reg > rho_rough * 1.01
                    results_exp2.reg_fail_delta(p_idx, atk_idx, trial) = true;
                end
                % 判定2: hinf_norm > (1/delta) * 1.01
                if hinf_reg > (1/delta_reg) * 1.01
                    results_exp2.reg_fail_hinf(p_idx, atk_idx, trial) = true;
                end
                
                fprintf('reg: delta=%.3e, hinf=%.3e\n', delta_reg, hinf_reg);
            catch ME
                fprintf('reg: INFEASIBLE\n');
                results_exp2.reg_delta(p_idx, atk_idx, trial) = NaN;
                results_exp2.reg_hinf(p_idx, atk_idx, trial) = NaN;
                results_exp2.reg_infeasible(p_idx, atk_idx, trial) = true;
                results_exp2.reg_fail_delta(p_idx, atk_idx, trial) = true;
                results_exp2.reg_fail_hinf(p_idx, atk_idx, trial) = true;
            end
        end
    end
end

% 攻撃成功率の計算
results_exp2.ori_success_rate_delta = squeeze(mean(results_exp2.ori_fail_delta, 3)) * 100;
results_exp2.ori_success_rate_hinf = squeeze(mean(results_exp2.ori_fail_hinf, 3)) * 100;
results_exp2.reg_success_rate_delta = squeeze(mean(results_exp2.reg_fail_delta, 3)) * 100;
results_exp2.reg_success_rate_hinf = squeeze(mean(results_exp2.reg_fail_hinf, 3)) * 100;

fprintf('\n=== 実験2 結果サマリ ===\n');
fprintf('攻撃成功率 (%%)\n');
fprintf('\n[original_thesis - 判定1: delta > rho_rough]\n');
fprintf('param\\attack\t');
for atk_idx = 1:length(attack_sizes)
    fprintf('%.0e\t\t', attack_sizes(atk_idx));
end
fprintf('\n');
for p_idx = 1:length(param_list)
    fprintf('%.0e\t\t', param_list(p_idx));
    for atk_idx = 1:length(attack_sizes)
        fprintf('%.1f%%\t\t', results_exp2.ori_success_rate_delta(p_idx, atk_idx));
    end
    fprintf('\n');
end

fprintf('\n[original_thesis - 判定2: hinf > 1/delta]\n');
fprintf('param\\attack\t');
for atk_idx = 1:length(attack_sizes)
    fprintf('%.0e\t\t', attack_sizes(atk_idx));
end
fprintf('\n');
for p_idx = 1:length(param_list)
    fprintf('%.0e\t\t', param_list(p_idx));
    for atk_idx = 1:length(attack_sizes)
        fprintf('%.1f%%\t\t', results_exp2.ori_success_rate_hinf(p_idx, atk_idx));
    end
    fprintf('\n');
end

fprintf('\n[regularization_sdp - 判定1: delta > rho_rough]\n');
fprintf('param\\attack\t');
for atk_idx = 1:length(attack_sizes)
    fprintf('%.0e\t\t', attack_sizes(atk_idx));
end
fprintf('\n');
for p_idx = 1:length(param_list)
    fprintf('%.0e\t\t', param_list(p_idx));
    for atk_idx = 1:length(attack_sizes)
        fprintf('%.1f%%\t\t', results_exp2.reg_success_rate_delta(p_idx, atk_idx));
    end
    fprintf('\n');
end

fprintf('\n[regularization_sdp - 判定2: hinf > 1/delta]\n');
fprintf('param\\attack\t');
for atk_idx = 1:length(attack_sizes)
    fprintf('%.0e\t\t', attack_sizes(atk_idx));
end
fprintf('\n');
for p_idx = 1:length(param_list)
    fprintf('%.0e\t\t', param_list(p_idx));
    for atk_idx = 1:length(attack_sizes)
        fprintf('%.1f%%\t\t', results_exp2.reg_success_rate_hinf(p_idx, atk_idx));
    end
    fprintf('\n');
end

%% ========== 結果保存 ==========
save('experiment_results.mat', 'results_exp1', 'results_exp2');
fprintf('\n結果を experiment_results.mat に保存しました。\n');

%% ========== グラフ描画 ==========
fprintf('\nグラフを描画中...\n');

% 実験1のグラフ
plot_experiment1(results_exp1);

% 実験2のグラフ
plot_experiment2(results_exp2);

fprintf('完了！\n');


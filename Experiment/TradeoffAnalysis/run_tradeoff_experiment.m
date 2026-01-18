% run_tradeoff_experiment.m
% 保守性（性能劣化）vs ロバスト性（攻撃成功率）のトレードオフ分析
%
% 散布図で2手法を比較
% - 横軸: 攻撃成功率（hinf判定）[%]
% - 縦軸: delta（大きいほど良い）
% - 色分け: original_thesis（青）、regularization_sdp（赤）
%
clear; clc; close all;

%% ========== パラメータ設定 ==========
% システム設定（cfg.Systemから取得）
A = cfg.System.A;
B = cfg.System.B;
[n, m] = cfg.System.getDimensions();
T = cfg.System.getSampleCount();

% パラメータ範囲
param_list = [1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1];

% 攻撃設定
attack_eps = 5e-3;                % 攻撃の大きさ（固定）
gamma_attack = 1e3;  % IDGSM攻撃時のgamma
gamma_reg = 1e3;     % regularization_sdpのgamma（固定）

% 試行回数
num_trials = 10;

% rho_roughの計算（cfg.Systemのヘルパー関数を使用）
rho_rough = cfg.System.calcRhoRough(A);
fprintf('rho_rough (理論的下界): %.6e\n', rho_rough);
fprintf('攻撃の大きさ: %.4f\n', attack_eps);
fprintf('パラメータ: ');
fprintf('%.0e ', param_list);
fprintf('\n');

%% ========== 結果格納 ==========
% 実験1: 保守性（攻撃なし）
results_exp1 = struct();
results_exp1.delta_ori = zeros(length(param_list), num_trials);
results_exp1.delta_reg = zeros(length(param_list), num_trials);
results_exp1.infeasible_ori = false(length(param_list), num_trials);
results_exp1.infeasible_reg = false(length(param_list), num_trials);

% 実験2: ロバスト性（攻撃あり）
results_exp2 = struct();
results_exp2.ori_fail_hinf = false(length(param_list), num_trials);
results_exp2.reg_fail_hinf = false(length(param_list), num_trials);
results_exp2.ori_infeasible = false(length(param_list), num_trials);
results_exp2.reg_infeasible = false(length(param_list), num_trials);

%% ========== 実験1: 保守性の評価（攻撃なし） ==========
fprintf('\n========== 実験1: 保守性の評価 ==========\n');

for trial = 1:num_trials
    fprintf('\n--- Trial %d/%d ---\n', trial, num_trials);
    
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
        
        data = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
        
        % --- original_thesis ---
        try
            [sol_ori, ~, ~, ~, ~] = baseline.solve_sdp(data);
            results_exp1.delta_ori(p_idx, trial) = sol_ori.delta;
            fprintf('ori=%.4e, ', sol_ori.delta);
        catch
            fprintf('ori=INFEASIBLE, ');
            results_exp1.delta_ori(p_idx, trial) = NaN;
            results_exp1.infeasible_ori(p_idx, trial) = true;
        end
        
        % --- regularization_sdp (with robust) ---
        try
            [~, ~, delta_reg, ~, ~] = proposed.solve_sdp_with_robust(data, param, gamma_reg);
            results_exp1.delta_reg(p_idx, trial) = delta_reg;
            fprintf('reg=%.4e\n', delta_reg);
        catch
            fprintf('reg=INFEASIBLE\n');
            results_exp1.delta_reg(p_idx, trial) = NaN;
            results_exp1.infeasible_reg(p_idx, trial) = true;
        end
    end
end

% 平均計算
results_exp1.delta_ori_mean = nanmean(results_exp1.delta_ori, 2);
results_exp1.delta_reg_mean = nanmean(results_exp1.delta_reg, 2);

%% ========== 実験2: ロバスト性の評価（攻撃あり） ==========
fprintf('\n========== 実験2: ロバスト性の評価 ==========\n');

% 基本のPhi行列（攻撃計算用）
Phi11_base = 1e-1 * eye(n);
Phi12 = zeros(n, T);
Phi22 = -eye(T);

for trial = 1:num_trials
    fprintf('\n===== Trial %d/%d =====\n', trial, num_trials);
    
    % データ生成（攻撃前）
    V = make_inputU(m, T);
    [X_ori, Z_ori] = datasim.simulate_openloop(A, B, V);
    U_ori = V;
    
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
        attack_success = true;
    catch ME
        fprintf('  [攻撃失敗] %s\n', ME.message);
        attack_success = false;
    end
    
    if ~attack_success
        % 攻撃失敗時は全パラメータをNaN
        for p_idx = 1:length(param_list)
            results_exp2.ori_fail_hinf(p_idx, trial) = NaN;
            results_exp2.reg_fail_hinf(p_idx, trial) = NaN;
        end
        continue;
    end
    
    % 各パラメータでSDP評価
    for p_idx = 1:length(param_list)
        param = param_list(p_idx);
        fprintf('  param = %.0e: ', param);
        
        Phi11 = param * eye(n);
        data_attacked = datasim.SystemData(A, B, X_atk, Z_atk, U_atk, Phi11, Phi12, Phi22);
        
        % --- original_thesis ---
        try
            [sol_ori, K_ori, ~, ~, ~] = baseline.solve_sdp(data_attacked);
            delta_ori = sol_ori.delta;
            hinf_ori = helper.hinfnorm_AK(A, B, K_ori);
            
            % 判定: hinf > 1/delta
            if hinf_ori > (1/delta_ori) * 1.01
                results_exp2.ori_fail_hinf(p_idx, trial) = true;
            end
            fprintf('ori: delta=%.3e, hinf=%.3e, ', delta_ori, hinf_ori);
        catch
            fprintf('ori: INFEASIBLE, ');
            results_exp2.ori_infeasible(p_idx, trial) = true;
            results_exp2.ori_fail_hinf(p_idx, trial) = true;
        end
        
        % --- regularization_sdp (with robust) ---
        try
            [sol_reg, K_reg, delta_reg, ~, ~] = proposed.solve_sdp_with_robust(data_attacked, param, gamma_reg);
            hinf_reg = helper.hinfnorm_AK(A, B, K_reg);
            
            % 判定: hinf > 1/delta
            if hinf_reg > (1/delta_reg) * 1.01
                results_exp2.reg_fail_hinf(p_idx, trial) = true;
            end
            fprintf('reg: delta=%.3e, hinf=%.3e\n', delta_reg, hinf_reg);
        catch
            fprintf('reg: INFEASIBLE\n');
            results_exp2.reg_infeasible(p_idx, trial) = true;
            results_exp2.reg_fail_hinf(p_idx, trial) = true;
        end
    end
end

% 攻撃成功率の計算（NaNを除外）
results_exp2.ori_success_rate = zeros(length(param_list), 1);
results_exp2.reg_success_rate = zeros(length(param_list), 1);

for p_idx = 1:length(param_list)
    valid_ori = ~isnan(results_exp2.ori_fail_hinf(p_idx, :));
    valid_reg = ~isnan(results_exp2.reg_fail_hinf(p_idx, :));
    
    if sum(valid_ori) > 0
        results_exp2.ori_success_rate(p_idx) = sum(results_exp2.ori_fail_hinf(p_idx, valid_ori)) / sum(valid_ori) * 100;
    else
        results_exp2.ori_success_rate(p_idx) = NaN;
    end
    
    if sum(valid_reg) > 0
        results_exp2.reg_success_rate(p_idx) = sum(results_exp2.reg_fail_hinf(p_idx, valid_reg)) / sum(valid_reg) * 100;
    else
        results_exp2.reg_success_rate(p_idx) = NaN;
    end
end

%% ========== 結果サマリ ==========
fprintf('\n========== 結果サマリ ==========\n');
fprintf('param\t\tori_delta\treg_delta\tori_attack%%\treg_attack%%\n');
for p_idx = 1:length(param_list)
    fprintf('%.0e\t\t%.4e\t%.4e\t%.1f%%\t\t%.1f%%\n', ...
        param_list(p_idx), ...
        results_exp1.delta_ori_mean(p_idx), ...
        results_exp1.delta_reg_mean(p_idx), ...
        results_exp2.ori_success_rate(p_idx), ...
        results_exp2.reg_success_rate(p_idx));
end

%% ========== 結果保存 ==========
% ファイル名に攻撃サイズを含める
filename_base = sprintf('tradeoff_eps%.0e', attack_eps);
save([filename_base '_results.mat'], 'results_exp1', 'results_exp2', 'param_list', 'rho_rough', 'attack_eps');
fprintf('\n結果を %s_results.mat に保存しました。\n', filename_base);

%% ========== 散布図描画 ==========
fprintf('\n散布図を描画中...\n');

figure('Name', 'Tradeoff Analysis', 'Position', [100, 100, 800, 600]);
hold on;

% original_thesis（青）
scatter(results_exp2.ori_success_rate, results_exp1.delta_ori_mean, ...
    150, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'original\_thesis');

% regularization_sdp（赤）
scatter(results_exp2.reg_success_rate, results_exp1.delta_reg_mean, ...
    150, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'regularization\_sdp');

% パラメータラベルを追加
for p_idx = 1:length(param_list)
    % original_thesis
    if ~isnan(results_exp2.ori_success_rate(p_idx)) && ~isnan(results_exp1.delta_ori_mean(p_idx))
        text(results_exp2.ori_success_rate(p_idx) + 2, results_exp1.delta_ori_mean(p_idx), ...
            sprintf('%.0e', param_list(p_idx)), 'Color', 'b', 'FontSize', 9);
    end
    % regularization_sdp
    if ~isnan(results_exp2.reg_success_rate(p_idx)) && ~isnan(results_exp1.delta_reg_mean(p_idx))
        text(results_exp2.reg_success_rate(p_idx) + 2, results_exp1.delta_reg_mean(p_idx), ...
            sprintf('%.0e', param_list(p_idx)), 'Color', 'r', 'FontSize', 9);
    end
end

% rho_roughの参照線
yline(rho_rough, '--', 'Color', [0.3, 0.7, 0.3], 'LineWidth', 2, ...
    'Label', sprintf('rho\\_rough = %.3e', rho_rough), 'LabelHorizontalAlignment', 'left');

hold off;

xlabel('攻撃成功率 [%]', 'FontSize', 14);
ylabel('delta（大きいほど良い）', 'FontSize', 14);
title(sprintf('保守性 vs ロバスト性のトレードオフ（攻撃サイズ: %.2f）', attack_eps), 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12);
xlim([-5, 105]);
grid on;
set(gca, 'FontSize', 12);

% 保存（ファイル名に攻撃サイズを含める）
saveas(gcf, [filename_base '_scatter.png']);
saveas(gcf, [filename_base '_scatter.fig']);
fprintf('散布図を %s_scatter.png に保存しました。\n', filename_base);

fprintf('完了！\n');


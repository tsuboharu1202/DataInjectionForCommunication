% demo_idgsm_normalization_effect.m
% IDGSMの正規化（normalization）の効果を評価
%
% 目的：
% 1. 10個のランダムシステムを生成（Aは不安定であることを確認）
% 2. 各システムに対して、IDGSM（正規化あり/なし）を実行
%    - 正規化あり: 勾配を要素単位の最大値で正規化（自然勾配法のような効果）
%    - 正規化なし: 勾配をそのまま使用
% 3. IDGSMの反復中にSDPがinfeasibleになったら、そのシステムを破棄して別のシステムを生成
% 4. 結果を比較して、正規化の効果を評価（伸び率と箱ひげ図）

clear; clc; close all;

fprintf('=== IDGSM正規化の効果評価テスト ===\n\n');

% ============================================
% パラメータ設定
% ============================================
n_systems = 10;  % システム数
[n, m] = cfg.System.getDimensions();
T = 20;  % サンプル数

% IDGSMパラメータ
eps_att = 1e-3;  % 攻撃強度
direction = 'negative';  % deltaを小さくする方向
gamma = 1e3;  % 正則化パラメータ

% Phi設定（cfg.Systemから取得）
Phi = cfg.System.getDefaultPhi(n, T);
Phi11 = Phi.Phi11;
Phi12 = Phi.Phi12;
Phi22 = Phi.Phi22;

% 結果を保存する配列
results = struct();
results.delta_init = [];
results.delta_final_with_norm = [];
results.delta_final_without_norm = [];
results.growth_rate_with_norm = [];      % 伸び率 = (delta_final - delta_init) / delta_init
results.growth_rate_without_norm = [];
results.success_with_norm = [];
results.success_without_norm = [];
results.system_idx = [];

% ============================================
% 各システムでテスト
% ============================================
fprintf('システム数: %d\n', n_systems);
fprintf('システムサイズ: n=%d, m=%d, T=%d\n\n', n, m, T);

sys_count = 0;
max_attempts_per_system = 10;  % システムあたりの最大試行回数

while sys_count < n_systems
    attempt = 0;
    system_valid = false;
    
    while ~system_valid && attempt < max_attempts_per_system
        attempt = attempt + 1;
        sys_count = sys_count + 1;
        
        fprintf('=== システム %d/%d (試行 %d) ===\n', sys_count, n_systems, attempt);
        
        % システム生成（不安定なシステムを生成）
        [A, B] = datasim.make_lti(n, m);
        
        % Aの固有値を確認
        A_eig = eig(A);
        max_real = max(real(A_eig));
        
        % 不安定でない場合は不安定化（実部が正になるようにシフト）
        if max_real < 0
            fprintf('  Aが安定（実部最大=%.4f）→不安定化します\n', max_real);
            A = A - (max_real - 0.1) * eye(size(A));  % 実部を正にする
            A_eig = eig(A);
            max_real = max(real(A_eig));
        end
        
        % 不安定であることを確認（実部が正であることを確認）
        if ~any(real(A_eig) > 0)
            fprintf('  システム %d: 不安定化に失敗しました。スキップします。\n', sys_count);
            continue;
        end
        
        fprintf('  Aの固有値（実部最大）: %.4f (不安定)\n', max_real);
        
        % データ生成
        V = make_inputU(m, T);
        [X, Z, U] = datasim.simulate_openloop_stable(A, B, V);
        
        data = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
        
        % 初期SDPを解く
        fprintf('  初期SDPを解く...\n');
        try
            [sol_init, ~, ~, ~, diagnostics_init] = proposed.solve_sdp(data, gamma);
            if diagnostics_init.problem ~= 0
                fprintf('    初期SDP解決失敗 (status=%d)。システムを破棄します。\n', diagnostics_init.problem);
                continue;
            end
            delta_init = sol_init.delta;
            fprintf('    初期delta = %.6e\n', delta_init);
        catch ME
            fprintf('    初期SDP解決エラー: %s。システムを破棄します。\n', ME.message);
            continue;
        end
        
        % IDGSM（正規化あり）を実行
        fprintf('  IDGSM（正規化あり）を実行...\n');
        opts_with_norm = struct();
        opts_with_norm.gamma = gamma;
        opts_with_norm.epsilon = eps_att;
        opts_with_norm.direction = direction;
        opts_with_norm.normalize_grad = true;
        
        try
            [X_idgsm_norm, Z_idgsm_norm, U_idgsm_norm, history_norm] = algorithms.idgsm_delta(data, opts_with_norm);
            
            % 最終SDPを解く
            data_final_norm = datasim.SystemData(A, B, X_idgsm_norm, Z_idgsm_norm, U_idgsm_norm, Phi11, Phi12, Phi22);
            [sol_final_norm, ~, ~, ~, diagnostics_final_norm] = proposed.solve_sdp(data_final_norm, gamma);
            
            if diagnostics_final_norm.problem ~= 0
                fprintf('    最終SDP解決失敗 (status=%d)。システムを破棄します。\n', diagnostics_final_norm.problem);
                continue;
            end
            
            delta_final_norm = sol_final_norm.delta;
            growth_rate_norm = (delta_final_norm - delta_init) / delta_init;
            results.success_with_norm(end+1) = true;
            fprintf('    最終delta = %.6e, 伸び率 = %.4f%%\n', delta_final_norm, 100 * growth_rate_norm);
            
        catch ME
            fprintf('    IDGSM（正規化あり）エラー: %s。システムを破棄します。\n', ME.message);
            continue;
        end
        
        % IDGSM（正規化なし）を実行
        fprintf('  IDGSM（正規化なし）を実行...\n');
        opts_without_norm = struct();
        opts_without_norm.gamma = gamma;
        opts_without_norm.epsilon = eps_att;
        opts_without_norm.direction = direction;
        opts_without_norm.normalize_grad = false;
        
        try
            [X_idgsm_no_norm, Z_idgsm_no_norm, U_idgsm_no_norm, history_no_norm] = algorithms.idgsm_delta(data, opts_without_norm);
            
            % 最終SDPを解く
            data_final_no_norm = datasim.SystemData(A, B, X_idgsm_no_norm, Z_idgsm_no_norm, U_idgsm_no_norm, Phi11, Phi12, Phi22);
            [sol_final_no_norm, ~, ~, ~, diagnostics_final_no_norm] = proposed.solve_sdp(data_final_no_norm, gamma);
            
            if diagnostics_final_no_norm.problem ~= 0
                fprintf('    最終SDP解決失敗 (status=%d)。システムを破棄します。\n', diagnostics_final_no_norm.problem);
                continue;
            end
            
            delta_final_no_norm = sol_final_no_norm.delta;
            growth_rate_no_norm = (delta_final_no_norm - delta_init) / delta_init;
            results.success_without_norm(end+1) = true;
            fprintf('    最終delta = %.6e, 伸び率 = %.4f%%\n', delta_final_no_norm, 100 * growth_rate_no_norm);
            
        catch ME
            fprintf('    IDGSM（正規化なし）エラー: %s。システムを破棄します。\n', ME.message);
            continue;
        end
        
        % 結果を保存
        results.delta_init(end+1) = delta_init;
        results.delta_final_with_norm(end+1) = delta_final_norm;
        results.delta_final_without_norm(end+1) = delta_final_no_norm;
        results.growth_rate_with_norm(end+1) = growth_rate_norm;
        results.growth_rate_without_norm(end+1) = growth_rate_no_norm;
        results.system_idx(end+1) = sys_count;
        
        system_valid = true;
        fprintf('  システム %d: 成功\n\n', sys_count);
    end
    
    if ~system_valid
        fprintf('  システム %d: 最大試行回数に達しました。スキップします。\n\n', sys_count);
    end
end

% ============================================
% 結果の統計
% ============================================
fprintf('\n=== 結果統計 ===\n');

n_valid = length(results.system_idx);

if n_valid == 0
    fprintf('有効な結果がありません。\n');
    return;
end

fprintf('有効なシステム数: %d\n\n', n_valid);

% 基本統計
fprintf('【正規化あり】\n');
fprintf('  伸び率平均: %.4f%%\n', 100 * mean(results.growth_rate_with_norm));
fprintf('  伸び率標準偏差: %.4f%%\n', 100 * std(results.growth_rate_with_norm));
fprintf('  伸び率最小: %.4f%%\n', 100 * min(results.growth_rate_with_norm));
fprintf('  伸び率最大: %.4f%%\n', 100 * max(results.growth_rate_with_norm));
fprintf('  成功率: %.1f%% (%d/%d)\n', ...
    100 * sum(results.success_with_norm) / n_valid, ...
    sum(results.success_with_norm), n_valid);

fprintf('\n【正規化なし】\n');
fprintf('  伸び率平均: %.4f%%\n', 100 * mean(results.growth_rate_without_norm));
fprintf('  伸び率標準偏差: %.4f%%\n', 100 * std(results.growth_rate_without_norm));
fprintf('  伸び率最小: %.4f%%\n', 100 * min(results.growth_rate_without_norm));
fprintf('  伸び率最大: %.4f%%\n', 100 * max(results.growth_rate_without_norm));
fprintf('  成功率: %.1f%% (%d/%d)\n', ...
    100 * sum(results.success_without_norm) / n_valid, ...
    sum(results.success_without_norm), n_valid);

% 比較
fprintf('\n【比較】\n');
growth_rate_diff = results.growth_rate_with_norm - results.growth_rate_without_norm;
fprintf('  伸び率差の平均（正規化あり - 正規化なし）: %.4f%%\n', 100 * mean(growth_rate_diff));
fprintf('  伸び率差の標準偏差: %.4f%%\n', 100 * std(growth_rate_diff));

% ============================================
% 結果をプロット
% ============================================
fprintf('\n=== 結果をプロット ===\n');

figure('Name', 'IDGSM正規化の効果', 'Position', [100, 100, 1200, 800]);

% 1. 伸び率の箱ひげ図
subplot(2, 2, 1);
boxplot([results.growth_rate_with_norm' * 100, results.growth_rate_without_norm' * 100], ...
    'Labels', {'正規化あり', '正規化なし'});
ylabel('伸び率 (%)');
title('伸び率の比較（箱ひげ図）');
grid on;

% 2. 伸び率の散布図
subplot(2, 2, 2);
plot(results.growth_rate_without_norm * 100, results.growth_rate_with_norm * 100, ...
    'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
plot([min([results.growth_rate_without_norm, results.growth_rate_with_norm]) * 100, ...
      max([results.growth_rate_without_norm, results.growth_rate_with_norm]) * 100], ...
     [min([results.growth_rate_without_norm, results.growth_rate_with_norm]) * 100, ...
      max([results.growth_rate_without_norm, results.growth_rate_with_norm]) * 100], ...
     'r--', 'LineWidth', 1);
xlabel('伸び率 (%) - 正規化なし');
ylabel('伸び率 (%) - 正規化あり');
title('伸び率の散布図');
legend('データ', 'y=x', 'Location', 'best');
grid on;

% 3. 伸び率差のヒストグラム
subplot(2, 2, 3);
histogram(growth_rate_diff * 100, 10);
xlabel('伸び率差 (%) (正規化あり - 正規化なし)');
ylabel('頻度');
title('伸び率差の分布');
hold on;
plot([0, 0], ylim, 'r--', 'LineWidth', 1);
grid on;

% 4. システムごとの伸び率比較
subplot(2, 2, 4);
plot(1:n_valid, results.growth_rate_with_norm * 100, 'b-o', ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '正規化あり');
hold on;
plot(1:n_valid, results.growth_rate_without_norm * 100, 'r--s', ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', '正規化なし');
xlabel('システム番号');
ylabel('伸び率 (%)');
title('システムごとの伸び率比較');
legend('Location', 'best');
grid on;

sgtitle(sprintf('IDGSM正規化の効果評価 (有効システム数: %d)', n_valid), ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('プロット完了\n');
fprintf('\n=== テスト完了 ===\n');


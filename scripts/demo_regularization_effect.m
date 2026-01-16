% demo_regularization_effect.m
% 正則化の有無がSDPの精度に与える影響を評価
%
% 目的：
% 1. 10個のランダムシステムを生成（Aは不安定でないことを確認）
% 2. 各システムに対して、正則化あり/なしでSDPを解く
% 3. 結果を比較して、正則化の効果を評価

clear; clc; close all;

fprintf('=== 正則化の効果評価テスト ===\n\n');

% ============================================
% パラメータ設定
% ============================================
n_systems = 10;  % システム数
[n, m, T] = deal(3, 1, cfg.Const.SAMPLE_COUNT);

% 正則化パラメータ
gamma_with_reg = 1e3;      % 正則化あり
gamma_without_reg = 0;    % 正則化なし（または非常に小さい値）

% Phi設定
Phi11 = 1e-1 * eye(n);
Phi12 = zeros(n, T);
Phi22 = -eye(T);

% 結果を保存する配列
results = struct();
results.delta_with_reg = [];
results.delta_without_reg = [];
results.success_with_reg = [];
results.success_without_reg = [];
results.A_eigenvalues = {};
results.system_idx = [];

% ============================================
% 各システムでテスト
% ============================================
fprintf('システム数: %d\n', n_systems);
fprintf('システムサイズ: n=%d, m=%d, T=%d\n\n', n, m, T);

for sys_idx = 1:n_systems
    fprintf('=== システム %d/%d ===\n', sys_idx, n_systems);
    
    % システム生成（不安定でないことを確認）
    [A, B] = datasim.make_lti(n, m);
    
    % Aの固有値を確認（不安定でないことを確認）
    A_eig = eig(A);
    max_real = max(real(A_eig));
    
    % 不安定な場合は安定化
    if max_real >= 0
        fprintf('  Aが不安定（実部最大=%.4f）→安定化します\n', max_real);
        A = A - (max_real + 0.1) * eye(size(A));
        A_eig = eig(A);
        max_real = max(real(A_eig));
    end
    
    % 安定性を確認（実部が負であることを確認）
    if ~all(real(A_eig) < 0)
        warning('システム %d: 安定化に失敗しました。スキップします。', sys_idx);
        continue;
    end
    
    results.A_eigenvalues{sys_idx} = A_eig;
    fprintf('  Aの固有値（実部最大）: %.4f (安定)\n', max_real);
    
    % データ生成
    V = make_inputU(m);
    [X, Z, U] = datasim.simulate_openloop_stable(A, B, V);
    
    data = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
    
    % 正則化ありでSDPを解く
    fprintf('  正則化あり (gamma=%.2e) でSDPを解く...\n', gamma_with_reg);
    try
        [sol_with_reg, ~, ~, ~, diagnostics_with_reg] = regularization_sdp.solve_sdp(data, gamma_with_reg);
        if diagnostics_with_reg.problem == 0
            delta_with_reg = sol_with_reg.delta;
            results.delta_with_reg(end+1) = delta_with_reg;
            results.success_with_reg(end+1) = true;
            fprintf('    delta = %.6e (成功)\n', delta_with_reg);
        else
            results.delta_with_reg(end+1) = NaN;
            results.success_with_reg(end+1) = false;
            fprintf('    SDP解決失敗 (status=%d)\n', diagnostics_with_reg.problem);
        end
    catch ME
        results.delta_with_reg(end+1) = NaN;
        results.success_with_reg(end+1) = false;
        fprintf('    エラー: %s\n', ME.message);
    end
    
    % 正則化なしでSDPを解く（gamma=0または非常に小さい値）
    fprintf('  正則化なし (gamma=%.2e) でSDPを解く...\n', gamma_without_reg);
    try
        [sol_without_reg, ~, ~, ~, diagnostics_without_reg] = regularization_sdp.solve_sdp(data, gamma_without_reg);
        if diagnostics_without_reg.problem == 0
            delta_without_reg = sol_without_reg.delta;
            results.delta_without_reg(end+1) = delta_without_reg;
            results.success_without_reg(end+1) = true;
            fprintf('    delta = %.6e (成功)\n', delta_without_reg);
        else
            results.delta_without_reg(end+1) = NaN;
            results.success_without_reg(end+1) = false;
            fprintf('    SDP解決失敗 (status=%d)\n', diagnostics_without_reg.problem);
        end
    catch ME
        results.delta_without_reg(end+1) = NaN;
        results.success_without_reg(end+1) = false;
        fprintf('    エラー: %s\n', ME.message);
    end
    
    results.system_idx(end+1) = sys_idx;
    fprintf('\n');
end

% ============================================
% 結果の統計
% ============================================
fprintf('\n=== 結果統計 ===\n');

% 成功したシステムのみで統計を計算
valid_idx = results.success_with_reg & results.success_without_reg;
n_valid = sum(valid_idx);

if n_valid == 0
    fprintf('有効な結果がありません。\n');
    return;
end

fprintf('有効なシステム数: %d/%d\n\n', n_valid, length(results.system_idx));

delta_with_reg_valid = results.delta_with_reg(valid_idx);
delta_without_reg_valid = results.delta_without_reg(valid_idx);

% 基本統計
fprintf('【正則化あり】\n');
fprintf('  delta平均: %.6e\n', mean(delta_with_reg_valid));
fprintf('  delta標準偏差: %.6e\n', std(delta_with_reg_valid));
fprintf('  delta最小: %.6e\n', min(delta_with_reg_valid));
fprintf('  delta最大: %.6e\n', max(delta_with_reg_valid));
fprintf('  成功率: %.1f%% (%d/%d)\n', ...
    100 * sum(results.success_with_reg) / length(results.system_idx), ...
    sum(results.success_with_reg), length(results.system_idx));

fprintf('\n【正則化なし】\n');
fprintf('  delta平均: %.6e\n', mean(delta_without_reg_valid));
fprintf('  delta標準偏差: %.6e\n', std(delta_without_reg_valid));
fprintf('  delta最小: %.6e\n', min(delta_without_reg_valid));
fprintf('  delta最大: %.6e\n', max(delta_without_reg_valid));
fprintf('  成功率: %.1f%% (%d/%d)\n', ...
    100 * sum(results.success_without_reg) / length(results.system_idx), ...
    sum(results.success_without_reg), length(results.system_idx));

% 比較
fprintf('\n【比較】\n');
delta_diff = delta_with_reg_valid - delta_without_reg_valid;
fprintf('  delta差の平均（正則化あり - 正則化なし）: %.6e\n', mean(delta_diff));
fprintf('  delta差の標準偏差: %.6e\n', std(delta_diff));
fprintf('  delta比の平均（正則化あり / 正則化なし）: %.4f\n', mean(delta_with_reg_valid ./ delta_without_reg_valid));
fprintf('  delta比の標準偏差: %.4f\n', std(delta_with_reg_valid ./ delta_without_reg_valid));

% ============================================
% 結果をプロット
% ============================================
fprintf('\n=== 結果をプロット ===\n');

figure('Name', '正則化の効果', 'Position', [100, 100, 1200, 800]);

% 1. deltaの比較（箱ひげ図）
subplot(2, 2, 1);
boxplot([delta_with_reg_valid', delta_without_reg_valid'], ...
    'Labels', {'正則化あり', '正則化なし'});
ylabel('delta');
title('deltaの比較（箱ひげ図）');
grid on;

% 2. deltaの散布図
subplot(2, 2, 2);
plot(delta_without_reg_valid, delta_with_reg_valid, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
plot([min([delta_without_reg_valid, delta_with_reg_valid]), ...
    max([delta_without_reg_valid, delta_with_reg_valid])], ...
    [min([delta_without_reg_valid, delta_with_reg_valid]), ...
    max([delta_without_reg_valid, delta_with_reg_valid])], ...
    'r--', 'LineWidth', 1);
xlabel('delta (正則化なし)');
ylabel('delta (正則化あり)');
title('deltaの散布図');
legend('データ', 'y=x', 'Location', 'best');
grid on;

% 3. delta差のヒストグラム
subplot(2, 2, 3);
histogram(delta_diff, 10);
xlabel('delta差 (正則化あり - 正則化なし)');
ylabel('頻度');
title('delta差の分布');
grid on;

% 4. delta比のヒストグラム
subplot(2, 2, 4);
delta_ratio = delta_with_reg_valid ./ delta_without_reg_valid;
histogram(delta_ratio, 10);
xlabel('delta比 (正則化あり / 正則化なし)');
ylabel('頻度');
title('delta比の分布');
hold on;
plot([1, 1], ylim, 'r--', 'LineWidth', 1);
grid on;

sgtitle(sprintf('正則化の効果評価 (有効システム数: %d/%d)', n_valid, length(results.system_idx)), ...
    'FontSize', 14, 'FontWeight', 'bold');

fprintf('プロット完了\n');
fprintf('\n=== テスト完了 ===\n');


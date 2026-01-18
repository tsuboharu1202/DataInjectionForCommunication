% demo_attack.m
% Regularization版のDGSMとIDGSMの比較
%
% 目的：
% 1. Regularization版のSDPを解く
% 2. DGSM（Direct）とIDGSM（Iterative）の攻撃を実行
% 3. negative/positive両方向でdeltaの変化を比較

clear; clc; close all;

fprintf('=== DGSM vs IDGSM 比較テスト (Regularization版) ===\n\n');

% ============================================
% 1. モデル作成とSDPを解く
% ============================================
fprintf('1. モデル作成とSDPを解く...\n');

% システム設定（cfg.Systemから取得）
A = cfg.System.A;
B = cfg.System.B;
[n, m] = cfg.System.getDimensions();
T = 20;  % サンプル数

disp('A');disp(A);
disp('B');disp(B);
disp('eig(A)');disp(eig(A));

fprintf('  システムサイズ: n=%d, m=%d, T=%d\n', n, m, T);

% 入力とデータ取得（stable版）
V = make_inputU(m, T);
[X, Z, U] = datasim.simulate_openloop_stable(A, B, V);
fprintf('  データ生成完了（stable版）\n');

% Phi設定（cfg.Systemから取得）
Phi = cfg.System.getDefaultPhi(n, T);
Phi11 = Phi.Phi11;
Phi12 = Phi.Phi12;
Phi22 = Phi.Phi22;

% 正則化付きSDPを解く
data = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
gamma = 1e3;  % 正則化パラメータ
[sol_init, ~, ~, ~, ~] = proposed.solve_sdp(data, gamma);

delta_init = sol_init.delta;
fprintf('  初期delta = %.6f\n', delta_init);

% ============================================
% 2. 攻撃パラメータ設定
% ============================================
fprintf('\n2. 攻撃パラメータ設定...\n');

eps_att = 1e-3;  % 攻撃強度

fprintf('  攻撃強度: eps_att = %.2e\n', eps_att);
fprintf('  gamma = %.2e\n', gamma);

% ============================================
% 3. DGSM攻撃（negative方向）
% ============================================
fprintf('\n3. DGSM攻撃（negative方向）...\n');

opts_dgsm_neg = struct('gamma', gamma, 'epsilon', eps_att, 'direction', 'negative');
[X_dgsm_neg, Z_dgsm_neg, U_dgsm_neg] = algorithms.dgsm_delta(data, opts_dgsm_neg);
data_dgsm_neg = datasim.SystemData(A, B, X_dgsm_neg, Z_dgsm_neg, U_dgsm_neg, Phi11, Phi12, Phi22);
[sol_dgsm_neg, ~, ~, ~, ~] = proposed.solve_sdp(data_dgsm_neg, gamma);
delta_dgsm_neg = sol_dgsm_neg.delta;

fprintf('  DGSM (negative): delta = %.6f (変化: %.6e, %.2f%%)\n', ...
    delta_dgsm_neg, delta_dgsm_neg - delta_init, ...
    100 * (delta_dgsm_neg - delta_init) / delta_init);

% ============================================
% 4. DGSM攻撃（positive方向）
% ============================================
fprintf('\n4. DGSM攻撃（positive方向）...\n');

opts_dgsm_pos = struct('gamma', gamma, 'epsilon', eps_att, 'direction', 'positive');
[X_dgsm_pos, Z_dgsm_pos, U_dgsm_pos] = algorithms.dgsm_delta(data, opts_dgsm_pos);
data_dgsm_pos = datasim.SystemData(A, B, X_dgsm_pos, Z_dgsm_pos, U_dgsm_pos, Phi11, Phi12, Phi22);
[sol_dgsm_pos, ~, ~, ~, ~] = proposed.solve_sdp(data_dgsm_pos, gamma);
delta_dgsm_pos = sol_dgsm_pos.delta;

fprintf('  DGSM (positive): delta = %.6f (変化: %.6e, %.2f%%)\n', ...
    delta_dgsm_pos, delta_dgsm_pos - delta_init, ...
    100 * (delta_dgsm_pos - delta_init) / delta_init);

% ============================================
% 5. IDGSM攻撃（negative方向）- 正規化あり
% ============================================
fprintf('\n5. IDGSM攻撃（negative方向）- 正規化あり...\n');

opts_idgsm_neg = struct();
opts_idgsm_neg.gamma = gamma;
opts_idgsm_neg.epsilon = eps_att;
opts_idgsm_neg.direction = 'negative';
opts_idgsm_neg.save_history = true;
opts_idgsm_neg.save_grad_norms = true;
opts_idgsm_neg.normalize_grad = true;

[X_idgsm_neg, Z_idgsm_neg, U_idgsm_neg, history_idgsm_neg] = algorithms.idgsm_delta(data, opts_idgsm_neg);
data_idgsm_neg = datasim.SystemData(A, B, X_idgsm_neg, Z_idgsm_neg, U_idgsm_neg, Phi11, Phi12, Phi22);
[sol_idgsm_neg, ~, ~, ~, ~] = proposed.solve_sdp(data_idgsm_neg, gamma);
delta_idgsm_neg = sol_idgsm_neg.delta;

fprintf('  IDGSM (negative, 正規化あり): delta = %.6f (変化: %.6e, %.2f%%)\n', ...
    delta_idgsm_neg, delta_idgsm_neg - delta_init, ...
    100 * (delta_idgsm_neg - delta_init) / delta_init);
fprintf('  反復回数: %d\n', history_idgsm_neg.iter_count);

% ============================================
% 5b. IDGSM攻撃（negative方向）- 正規化なし
% ============================================
fprintf('\n5b. IDGSM攻撃（negative方向）- 正規化なし...\n');

opts_idgsm_neg_no_norm = struct();
opts_idgsm_neg_no_norm.gamma = gamma;
opts_idgsm_neg_no_norm.epsilon = eps_att;
opts_idgsm_neg_no_norm.direction = 'negative';
opts_idgsm_neg_no_norm.save_history = true;
opts_idgsm_neg_no_norm.save_grad_norms = true;
opts_idgsm_neg_no_norm.normalize_grad = false;

[X_idgsm_neg_no_norm, Z_idgsm_neg_no_norm, U_idgsm_neg_no_norm, history_idgsm_neg_no_norm] = algorithms.idgsm_delta(data, opts_idgsm_neg_no_norm);
data_idgsm_neg_no_norm = datasim.SystemData(A, B, X_idgsm_neg_no_norm, Z_idgsm_neg_no_norm, U_idgsm_neg_no_norm, Phi11, Phi12, Phi22);
[sol_idgsm_neg_no_norm, ~, ~, ~, ~] = proposed.solve_sdp(data_idgsm_neg_no_norm, gamma);
delta_idgsm_neg_no_norm = sol_idgsm_neg_no_norm.delta;

fprintf('  IDGSM (negative, 正規化なし): delta = %.6f (変化: %.6e, %.2f%%)\n', ...
    delta_idgsm_neg_no_norm, delta_idgsm_neg_no_norm - delta_init, ...
    100 * (delta_idgsm_neg_no_norm - delta_init) / delta_init);
fprintf('  反復回数: %d\n', history_idgsm_neg_no_norm.iter_count);

% ============================================
% 6. IDGSM攻撃（positive方向）
% ============================================
fprintf('\n6. IDGSM攻撃（positive方向）...\n');

opts_idgsm_pos = struct();
opts_idgsm_pos.gamma = gamma;
opts_idgsm_pos.epsilon = eps_att;
opts_idgsm_pos.direction = 'positive';
opts_idgsm_pos.save_history = true;
opts_idgsm_pos.save_grad_norms = true;
opts_idgsm_pos.normalize_grad = true;

[X_idgsm_pos, Z_idgsm_pos, U_idgsm_pos, history_idgsm_pos] = algorithms.idgsm_delta(data, opts_idgsm_pos);
data_idgsm_pos = datasim.SystemData(A, B, X_idgsm_pos, Z_idgsm_pos, U_idgsm_pos, Phi11, Phi12, Phi22);
[sol_idgsm_pos, ~, ~, ~, ~] = proposed.solve_sdp(data_idgsm_pos, gamma);
delta_idgsm_pos = sol_idgsm_pos.delta;

fprintf('  IDGSM (positive): delta = %.6f (変化: %.6e, %.2f%%)\n', ...
    delta_idgsm_pos, delta_idgsm_pos - delta_init, ...
    100 * (delta_idgsm_pos - delta_init) / delta_init);
fprintf('  反復回数: %d\n', history_idgsm_pos.iter_count);

% ============================================
% 7. 結果比較
% ============================================
fprintf('\n=== 結果比較 ===\n');
fprintf('初期delta: %.6f\n\n', delta_init);

fprintf('【Negative方向（deltaを小さくする）】\n');
fprintf('  DGSM:  delta = %.6f, 変化量 = %.6e (%.2f%%)\n', ...
    delta_dgsm_neg, delta_dgsm_neg - delta_init, ...
    100 * (delta_dgsm_neg - delta_init) / delta_init);
fprintf('  IDGSM (正規化あり): delta = %.6f, 変化量 = %.6e (%.2f%%)\n', ...
    delta_idgsm_neg, delta_idgsm_neg - delta_init, ...
    100 * (delta_idgsm_neg - delta_init) / delta_init);
fprintf('  IDGSM (正規化なし): delta = %.6f, 変化量 = %.6e (%.2f%%)\n', ...
    delta_idgsm_neg_no_norm, delta_idgsm_neg_no_norm - delta_init, ...
    100 * (delta_idgsm_neg_no_norm - delta_init) / delta_init);
fprintf('  差異: IDGSM(正規化あり) - DGSM = %.6e\n', delta_idgsm_neg - delta_dgsm_neg);
fprintf('  差異: IDGSM(正規化なし) - DGSM = %.6e\n', delta_idgsm_neg_no_norm - delta_dgsm_neg);
fprintf('  差異: IDGSM(正規化あり) - IDGSM(正規化なし) = %.6e\n\n', delta_idgsm_neg - delta_idgsm_neg_no_norm);

fprintf('【Positive方向（deltaを大きくする）】\n');
fprintf('  DGSM:  delta = %.6f, 変化量 = %.6e (%.2f%%)\n', ...
    delta_dgsm_pos, delta_dgsm_pos - delta_init, ...
    100 * (delta_dgsm_pos - delta_init) / delta_init);
fprintf('  IDGSM: delta = %.6f, 変化量 = %.6e (%.2f%%)\n', ...
    delta_idgsm_pos, delta_idgsm_pos - delta_init, ...
    100 * (delta_idgsm_pos - delta_init) / delta_init);
fprintf('  差異: IDGSM - DGSM = %.6e\n\n', delta_idgsm_pos - delta_dgsm_pos);

% ============================================
% 8. 勾配ノルムの変化をプロット（正規化あり/なしの比較）
% ============================================
fprintf('\n=== 勾配ノルムの変化をプロット（正規化あり/なしの比較） ===\n');

% Negative方向の勾配ノルムをプロット（正規化あり/なしの比較）
if isfield(history_idgsm_neg, 'grad_norms_X_Linf') && isfield(history_idgsm_neg_no_norm, 'grad_norms_X_Linf')
    fprintf('IDGSM (negative方向) - 正規化あり/なしの比較をプロット中...\n');
    
    iter_neg = 0:(history_idgsm_neg.iter_count-1);
    iter_neg_no_norm = 0:(history_idgsm_neg_no_norm.iter_count-1);
    
    figure('Name', 'IDGSM (negative) - 正規化あり/なしの比較', 'Position', [100, 100, 1400, 900]);
    
    % X勾配 - L∞ノルム
    subplot(3, 3, 1);
    semilogy(iter_neg, history_idgsm_neg.grad_norms_X_Linf, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化あり');
    hold on;
    semilogy(iter_neg_no_norm, history_idgsm_neg_no_norm.grad_norms_X_Linf, 'b--s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化なし');
    xlabel('反復回数');
    ylabel('L∞ノルム');
    title('X勾配のL∞ノルム');
    legend('Location', 'best');
    grid on;
    
    % X勾配 - Frobeniusノルム
    subplot(3, 3, 2);
    semilogy(iter_neg, history_idgsm_neg.grad_norms_X_Fro, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化あり');
    hold on;
    semilogy(iter_neg_no_norm, history_idgsm_neg_no_norm.grad_norms_X_Fro, 'b--s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化なし');
    xlabel('反復回数');
    ylabel('Frobeniusノルム');
    title('X勾配のFrobeniusノルム');
    legend('Location', 'best');
    grid on;
    
    % X勾配 - L∞/Frobenius比
    subplot(3, 3, 3);
    ratio_X_norm = history_idgsm_neg.grad_norms_X_Linf ./ history_idgsm_neg.grad_norms_X_Fro;
    ratio_X_no_norm = history_idgsm_neg_no_norm.grad_norms_X_Linf ./ history_idgsm_neg_no_norm.grad_norms_X_Fro;
    plot(iter_neg, ratio_X_norm, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化あり');
    hold on;
    plot(iter_neg_no_norm, ratio_X_no_norm, 'b--s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化なし');
    xlabel('反復回数');
    ylabel('L∞/Frobenius比');
    title('X勾配: L∞/Frobenius比（1に近い=1要素のみ大きい）');
    legend('Location', 'best');
    grid on;
    
    % Z勾配 - L∞ノルム
    subplot(3, 3, 4);
    semilogy(iter_neg, history_idgsm_neg.grad_norms_Z_Linf, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化あり');
    hold on;
    semilogy(iter_neg_no_norm, history_idgsm_neg_no_norm.grad_norms_Z_Linf, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化なし');
    xlabel('反復回数');
    ylabel('L∞ノルム');
    title('Z勾配のL∞ノルム');
    legend('Location', 'best');
    grid on;
    
    % Z勾配 - Frobeniusノルム
    subplot(3, 3, 5);
    semilogy(iter_neg, history_idgsm_neg.grad_norms_Z_Fro, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化あり');
    hold on;
    semilogy(iter_neg_no_norm, history_idgsm_neg_no_norm.grad_norms_Z_Fro, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化なし');
    xlabel('反復回数');
    ylabel('Frobeniusノルム');
    title('Z勾配のFrobeniusノルム');
    legend('Location', 'best');
    grid on;
    
    % Z勾配 - L∞/Frobenius比
    subplot(3, 3, 6);
    ratio_Z_norm = history_idgsm_neg.grad_norms_Z_Linf ./ history_idgsm_neg.grad_norms_Z_Fro;
    ratio_Z_no_norm = history_idgsm_neg_no_norm.grad_norms_Z_Linf ./ history_idgsm_neg_no_norm.grad_norms_Z_Fro;
    plot(iter_neg, ratio_Z_norm, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化あり');
    hold on;
    plot(iter_neg_no_norm, ratio_Z_no_norm, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化なし');
    xlabel('反復回数');
    ylabel('L∞/Frobenius比');
    title('Z勾配: L∞/Frobenius比（1に近い=1要素のみ大きい）');
    legend('Location', 'best');
    grid on;
    
    % U勾配 - L∞ノルム
    subplot(3, 3, 7);
    semilogy(iter_neg, history_idgsm_neg.grad_norms_U_Linf, 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化あり');
    hold on;
    semilogy(iter_neg_no_norm, history_idgsm_neg_no_norm.grad_norms_U_Linf, 'g--s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化なし');
    xlabel('反復回数');
    ylabel('L∞ノルム');
    title('U勾配のL∞ノルム');
    legend('Location', 'best');
    grid on;
    
    % U勾配 - Frobeniusノルム
    subplot(3, 3, 8);
    semilogy(iter_neg, history_idgsm_neg.grad_norms_U_Fro, 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化あり');
    hold on;
    semilogy(iter_neg_no_norm, history_idgsm_neg_no_norm.grad_norms_U_Fro, 'g--s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化なし');
    xlabel('反復回数');
    ylabel('Frobeniusノルム');
    title('U勾配のFrobeniusノルム');
    legend('Location', 'best');
    grid on;
    
    % U勾配 - L∞/Frobenius比
    subplot(3, 3, 9);
    ratio_U_norm = history_idgsm_neg.grad_norms_U_Linf ./ history_idgsm_neg.grad_norms_U_Fro;
    ratio_U_no_norm = history_idgsm_neg_no_norm.grad_norms_U_Linf ./ history_idgsm_neg_no_norm.grad_norms_U_Fro;
    plot(iter_neg, ratio_U_norm, 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化あり');
    hold on;
    plot(iter_neg_no_norm, ratio_U_no_norm, 'g--s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '正規化なし');
    xlabel('反復回数');
    ylabel('L∞/Frobenius比');
    title('U勾配: L∞/Frobenius比（1に近い=1要素のみ大きい）');
    legend('Location', 'best');
    grid on;
    
    sgtitle('IDGSM (negative方向) - 正規化あり/なしの比較', 'FontSize', 14, 'FontWeight', 'bold');
end

% 勾配の要素分布をプロット
if isfield(history_idgsm_neg, 'grad_max_ratio_X')
    fprintf('勾配の要素分布をプロット中...\n');
    
    figure('Name', '勾配の要素分布', 'Position', [200, 200, 1200, 600]);
    
    iter_neg = 0:(history_idgsm_neg.iter_count-1);
    
    % 最大要素の割合
    subplot(2, 3, 1);
    plot(iter_neg, history_idgsm_neg.grad_max_ratio_X, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlabel('反復回数');
    ylabel('最大要素の割合');
    title('X勾配: 最大要素/全要素和');
    grid on;
    ylim([0, 1]);
    
    subplot(2, 3, 2);
    plot(iter_neg, history_idgsm_neg.grad_max_ratio_Z, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlabel('反復回数');
    ylabel('最大要素の割合');
    title('Z勾配: 最大要素/全要素和');
    grid on;
    ylim([0, 1]);
    
    subplot(2, 3, 3);
    plot(iter_neg, history_idgsm_neg.grad_max_ratio_U, 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlabel('反復回数');
    ylabel('最大要素の割合');
    title('U勾配: 最大要素/全要素和');
    grid on;
    ylim([0, 1]);
    
    % 上位3要素の割合
    subplot(2, 3, 4);
    plot(iter_neg, history_idgsm_neg.grad_top3_ratio_X, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlabel('反復回数');
    ylabel('上位3要素の割合');
    title('X勾配: 上位3要素/全要素和');
    grid on;
    ylim([0, 1]);
    
    subplot(2, 3, 5);
    plot(iter_neg, history_idgsm_neg.grad_top3_ratio_Z, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlabel('反復回数');
    ylabel('上位3要素の割合');
    title('Z勾配: 上位3要素/全要素和');
    grid on;
    ylim([0, 1]);
    
    subplot(2, 3, 6);
    plot(iter_neg, history_idgsm_neg.grad_top3_ratio_U, 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlabel('反復回数');
    ylabel('上位3要素の割合');
    title('U勾配: 上位3要素/全要素和');
    grid on;
    ylim([0, 1]);
    
    sgtitle('勾配の要素分布（正規化あり）', 'FontSize', 14, 'FontWeight', 'bold');
end

fprintf('プロット完了\n');
fprintf('\n=== テスト完了 ===\n');

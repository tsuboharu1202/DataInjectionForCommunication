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

% システムサイズ
[n, m, T] = deal(3, 1, cfg.Const.SAMPLE_COUNT);

% システム定義（demo_implicit_regularization.mと同じ）
A = [-0.192, -0.936, -0.814;
    -0.918,  0.729, -0.724;
    -0.412,  0.735, -0.516];
B = [-0.554;
    0.735;
    0.528];

fprintf('  システムサイズ: n=%d, m=%d, T=%d\n', n, m, T);

% 入力とデータ取得（stable版）
V = make_inputU(m);
[X, Z, U] = datasim.simulate_openloop_stable(A, B, V);
fprintf('  データ生成完了（stable版）\n');

% Phi設定
Phi11 = 1e-1 * eye(n);
Phi12 = zeros(n, T);
Phi22 = -eye(T);

% 正則化付きSDPを解く
data = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
gamma = 1e3;
[sol_init, ~, ~, ~, ~] = regularization_sdp.solve_sdp(data, gamma);

delta_init = sol_init.delta;
fprintf('  初期delta = %.6f\n', delta_init);

% ============================================
% 2. 攻撃パラメータ設定
% ============================================
fprintf('\n2. 攻撃パラメータ設定...\n');

eps_att = cfg.Const.ATTACKER_UPPERLIMIT;
opts = struct('gamma', gamma);

fprintf('  攻撃強度: eps_att = %.2e\n', eps_att);
fprintf('  gamma = %.2e\n', gamma);

% ============================================
% 3. DGSM攻撃（negative方向）
% ============================================
fprintf('\n3. DGSM攻撃（negative方向）...\n');

[X_dgsm_neg, Z_dgsm_neg, U_dgsm_neg] = attack.dgsm_delta(data, eps_att, 'negative', opts);
data_dgsm_neg = datasim.SystemData(A, B, X_dgsm_neg, Z_dgsm_neg, U_dgsm_neg, Phi11, Phi12, Phi22);
[sol_dgsm_neg, ~, ~, ~, ~] = regularization_sdp.solve_sdp(data_dgsm_neg, gamma);
delta_dgsm_neg = sol_dgsm_neg.delta;

fprintf('  DGSM (negative): delta = %.6f (変化: %.6e, %.2f%%)\n', ...
    delta_dgsm_neg, delta_dgsm_neg - delta_init, ...
    100 * (delta_dgsm_neg - delta_init) / delta_init);

% ============================================
% 4. DGSM攻撃（positive方向）
% ============================================
fprintf('\n4. DGSM攻撃（positive方向）...\n');

[X_dgsm_pos, Z_dgsm_pos, U_dgsm_pos] = attack.dgsm_delta(data, eps_att, 'positive', opts);
data_dgsm_pos = datasim.SystemData(A, B, X_dgsm_pos, Z_dgsm_pos, U_dgsm_pos, Phi11, Phi12, Phi22);
[sol_dgsm_pos, ~, ~, ~, ~] = regularization_sdp.solve_sdp(data_dgsm_pos, gamma);
delta_dgsm_pos = sol_dgsm_pos.delta;

fprintf('  DGSM (positive): delta = %.6f (変化: %.6e, %.2f%%)\n', ...
    delta_dgsm_pos, delta_dgsm_pos - delta_init, ...
    100 * (delta_dgsm_pos - delta_init) / delta_init);

% ============================================
% 5. IDGSM攻撃（negative方向）
% ============================================
fprintf('\n5. IDGSM攻撃（negative方向）...\n');

[X_idgsm_neg, Z_idgsm_neg, U_idgsm_neg] = attack.idgsm_delta(data, false, [], [], eps_att, 'negative', opts);
data_idgsm_neg = datasim.SystemData(A, B, X_idgsm_neg, Z_idgsm_neg, U_idgsm_neg, Phi11, Phi12, Phi22);
[sol_idgsm_neg, ~, ~, ~, ~] = regularization_sdp.solve_sdp(data_idgsm_neg, gamma);
delta_idgsm_neg = sol_idgsm_neg.delta;

fprintf('  IDGSM (negative): delta = %.6f (変化: %.6e, %.2f%%)\n', ...
    delta_idgsm_neg, delta_idgsm_neg - delta_init, ...
    100 * (delta_idgsm_neg - delta_init) / delta_init);

% ============================================
% 6. IDGSM攻撃（positive方向）
% ============================================
fprintf('\n6. IDGSM攻撃（positive方向）...\n');

[X_idgsm_pos, Z_idgsm_pos, U_idgsm_pos] = attack.idgsm_delta(data, false, [], [], eps_att, 'positive', opts);
data_idgsm_pos = datasim.SystemData(A, B, X_idgsm_pos, Z_idgsm_pos, U_idgsm_pos, Phi11, Phi12, Phi22);
[sol_idgsm_pos, ~, ~, ~, ~] = regularization_sdp.solve_sdp(data_idgsm_pos, gamma);
delta_idgsm_pos = sol_idgsm_pos.delta;

fprintf('  IDGSM (positive): delta = %.6f (変化: %.6e, %.2f%%)\n', ...
    delta_idgsm_pos, delta_idgsm_pos - delta_init, ...
    100 * (delta_idgsm_pos - delta_init) / delta_init);

% ============================================
% 7. 結果比較
% ============================================
fprintf('\n=== 結果比較 ===\n');
fprintf('初期delta: %.6f\n\n', delta_init);

fprintf('【Negative方向（deltaを小さくする）】\n');
fprintf('  DGSM:  delta = %.6f, 変化量 = %.6e (%.2f%%)\n', ...
    delta_dgsm_neg, delta_dgsm_neg - delta_init, ...
    100 * (delta_dgsm_neg - delta_init) / delta_init);
fprintf('  IDGSM: delta = %.6f, 変化量 = %.6e (%.2f%%)\n', ...
    delta_idgsm_neg, delta_idgsm_neg - delta_init, ...
    100 * (delta_idgsm_neg - delta_init) / delta_init);
fprintf('  差異: IDGSM - DGSM = %.6e\n\n', delta_idgsm_neg - delta_dgsm_neg);

fprintf('【Positive方向（deltaを大きくする）】\n');
fprintf('  DGSM:  delta = %.6f, 変化量 = %.6e (%.2f%%)\n', ...
    delta_dgsm_pos, delta_dgsm_pos - delta_init, ...
    100 * (delta_dgsm_pos - delta_init) / delta_init);
fprintf('  IDGSM: delta = %.6f, 変化量 = %.6e (%.2f%%)\n', ...
    delta_idgsm_pos, delta_idgsm_pos - delta_init, ...
    100 * (delta_idgsm_pos - delta_init) / delta_init);
fprintf('  差異: IDGSM - DGSM = %.6e\n\n', delta_idgsm_pos - delta_dgsm_pos);

fprintf('=== テスト完了 ===\n');

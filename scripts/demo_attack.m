% U生成→データ生成→SDP→表示の最小デモ
clear; clc; close all;

% 1) 連続の種 or 再現性
% rng(1);

[n,m,T] = deal(5,3,cfg.Const.SAMPLE_COUNT);
[A,B] = datasim.make_lti(n,m);

disp('A');disp(A);
disp('B');disp(B);

% 3) 入力とデータ取得
V = make_inputU(m);
[X,Z,U] = datasim.simulate_openloop_stable(A,B,V);


W = Z - A*X - B*U;   % n×T

Phi11 = 1e-7*eye(n);     % 5%マージン
% Phi11 = cfg.Const.SAMPLE_COUNT * eye(n);
Phi12 = zeros(n,T);
Phi22 = -eye(T);


% 4) SDP (regularization版)
data = datasim.SystemData(A,B,X,Z,U,Phi11,Phi12,Phi22);
gamma = 1e3;
[sol, K, ~, ~, ~] = regularization_sdp.solve_sdp(data, gamma);
K_ori = K;
% Use full eig for robustness (small n), and compute spectral radius.
ev_ori = eig(A + B*K_ori);
rho_ori = max(abs(ev_ori));
delta_ori = sol.delta;
fprintf('元のrho: %e, delta: %e\n', rho_ori, delta_ori);


% ============================================
% 攻撃の実行
% ============================================
fprintf('\n=== 攻撃の実行 ===\n');

% 1. DIRECT_DGSM_DELTA (deltaを小さくする方向)
fprintf('\n1. DIRECT_DGSM_DELTA 攻撃（deltaを小さくする方向）...\n');
[X_dgsm, Z_dgsm, U_dgsm] = attack.execute_attack(data, cfg.AttackType.DIRECT_DGSM_DELTA);
sd_dgsm = datasim.SystemData(A,B,X_dgsm,Z_dgsm,U_dgsm,data.Phi11,data.Phi12,data.Phi22);
gamma = 1e3;
[sol_dgsm, ~, ~, ~, ~] = regularization_sdp.solve_sdp(sd_dgsm, gamma);
delta_dgsm = sol_dgsm.delta;
K_dgsm = sol_dgsm.K;
rho_dgsm = max(abs(eig(A+B*K_dgsm)));
fprintf('  delta_dgsm: %e (元: %e, 変化: %+.2e)\n', delta_dgsm, delta_ori, delta_dgsm - delta_ori);
fprintf('  rho_dgsm: %e (元: %e, 変化: %+.2e)\n', rho_dgsm, rho_ori, rho_dgsm - rho_ori);

% 2. IMPLICIT_IDGSM_DELTA_POSITIVE (deltaを大きくする方向)
fprintf('\n2. IMPLICIT_IDGSM_DELTA_POSITIVE 攻撃（deltaを大きくする方向）...\n');
[X_large, Z_large, U_large, history_large] = attack.execute_attack(data, cfg.AttackType.IMPLICIT_IDGSM_DELTA_POSITIVE, [], [], true);
sd_large = datasim.SystemData(A,B,X_large,Z_large,U_large,data.Phi11,data.Phi12,data.Phi22);
[sol_large, ~, ~, ~, ~] = regularization_sdp.solve_sdp(sd_large, gamma);
delta_large = sol_large.delta;
K_large = sol_large.K;
rho_large = max(abs(eig(A+B*K_large)));
fprintf('  delta_large: %e (元: %e, 変化: %+.2e)\n', delta_large, delta_ori, delta_large - delta_ori);
fprintf('  rho_large: %e (元: %e, 変化: %+.2e)\n', rho_large, rho_ori, rho_large - rho_ori);
if ~isempty(history_large) && isfield(history_large, 'delta_history')
    fprintf('  履歴: delta[0]=%e -> delta[%d]=%e (反復回数: %d)\n', ...
        history_large.delta_history(1), length(history_large.delta_history)-1, ...
        history_large.delta_history(end), history_large.iter_count);
end

% 3. IMPLICIT_IDGSM_DELTA_NEGATIVE (deltaを小さくする方向)
fprintf('\n3. IMPLICIT_IDGSM_DELTA_NEGATIVE 攻撃（deltaを小さくする方向）...\n');
[X_small, Z_small, U_small, history_small] = attack.execute_attack(data, cfg.AttackType.IMPLICIT_IDGSM_DELTA_NEGATIVE, [], [], true);
sd_small = datasim.SystemData(A,B,X_small,Z_small,U_small,data.Phi11,data.Phi12,data.Phi22);
[sol_small, ~, ~, ~, ~] = regularization_sdp.solve_sdp(sd_small, gamma);
delta_small = sol_small.delta;
K_small = sol_small.K;
rho_small = max(abs(eig(A+B*K_small)));
fprintf('  delta_small: %e (元: %e, 変化: %+.2e)\n', delta_small, delta_ori, delta_small - delta_ori);
fprintf('  rho_small: %e (元: %e, 変化: %+.2e)\n', rho_small, rho_ori, rho_small - rho_ori);
if ~isempty(history_small) && isfield(history_small, 'delta_history')
    fprintf('  履歴: delta[0]=%e -> delta[%d]=%e (反復回数: %d)\n', ...
        history_small.delta_history(1), length(history_small.delta_history)-1, ...
        history_small.delta_history(end), history_small.iter_count);
end

% ============================================
% 結果の比較
% ============================================
fprintf('\n=== 結果の比較 ===\n');
fprintf('元のdelta: %e, rho: %e\n', delta_ori, rho_ori);
fprintf('\n--- Deltaの変化 ---\n');
fprintf('DIRECT_DGSM_DELTA:        delta=%e (変化: %+.2e)\n', delta_dgsm, delta_dgsm - delta_ori);
fprintf('IMPLICIT_IDGSM_DELTA_POSITIVE:  delta=%e (変化: %+.2e)\n', delta_large, delta_large - delta_ori);
fprintf('IMPLICIT_IDGSM_DELTA_NEGATIVE:  delta=%e (変化: %+.2e)\n', delta_small, delta_small - delta_ori);
fprintf('\n--- Rhoの変化 ---\n');
fprintf('DIRECT_DGSM_DELTA:        rho=%e (変化: %+.2e)\n', rho_dgsm, rho_dgsm - rho_ori);
fprintf('IMPLICIT_IDGSM_DELTA_POSITIVE:  rho=%e (変化: %+.2e)\n', rho_large, rho_large - rho_ori);
fprintf('IMPLICIT_IDGSM_DELTA_NEGATIVE:  rho=%e (変化: %+.2e)\n', rho_small, rho_small - rho_ori);

% ============================================
% 可視化
% ============================================
fprintf('\n=== 可視化 ===\n');
fprintf('DIRECT_DGSM_DELTAの結果を表示します...\n');
visualize.plot_data(X, X_dgsm);
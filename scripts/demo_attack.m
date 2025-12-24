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


% 4) SDP
data = datasim.SystemData(A,B,X,Z,U,Phi11,Phi12,Phi22);
[sol, K, Y, L, diagnostics] = data.solve_sdp_on_data();
K_ori =  K;
% Use full eig for robustness (small n), and compute spectral radius.
ev_ori = eig(A + B*K_ori);
rho_ori = max(abs(ev_ori));
disp('rho_ori');disp(rho_ori);


% ============================================
% 攻撃の実行
% ============================================
fprintf('\n=== 攻撃の実行 ===\n');

% % 1. DIRECT_DGSM_EV
fprintf('\n1. DIRECT_DGSM_EV 攻撃...\n');
[X_dgsm, Z_dgsm, U_dgsm] = attack.execute_attack(data, cfg.AttackType.DIRECT_DGSM_EV);
sd_dgsm = datasim.SystemData(A,B,X_dgsm,Z_dgsm,U_dgsm,data.Phi11,data.Phi12,data.Phi22);
[~, K_dgsm, ~, ~, ~] = sd_dgsm.solve_sdp_on_data();
rho_dgsm = max(abs(eig(A+B*K_dgsm)));
fprintf('  rho_dgsm: %e (元: %e)\n', rho_dgsm, rho_ori);

% 2. IMPLICIT_IDGSM_RHO_LARGE (勾配と同じ方向、Deltaを大きくする)
fprintf('\n2. IMPLICIT_IDGSM_RHO_LARGE 攻撃（勾配と同じ方向）...\n');
[X_large, Z_large, U_large, history_large] = attack.execute_attack(data, cfg.AttackType.IMPLICIT_IDGSM_RHO_LARGE, [], 10, true);
sd_large = datasim.SystemData(A,B,X_large,Z_large,U_large,data.Phi11,data.Phi12,data.Phi22);
[~, K_large, ~, ~, ~] = sd_large.solve_sdp_on_data();
rho_large = max(abs(eig(A+B*K_large)));
fprintf('  rho_large: %e (元: %e)\n', rho_large, rho_ori);
if ~isempty(history_large) && isfield(history_large, 'rho')
    fprintf('  履歴: rho[0]=%e -> rho[%d]=%e\n', history_large.rho(1), length(history_large.rho)-1, history_large.rho(end));
end

% 3. IMPLICIT_IDGSM_RHO_SMALL (勾配と逆方向、Deltaを小さくする)
fprintf('\n3. IMPLICIT_IDGSM_RHO_SMALL 攻撃（勾配と逆方向）...\n');
[X_small, Z_small, U_small, history_small] = attack.execute_attack(data, cfg.AttackType.IMPLICIT_IDGSM_RHO_SMALL, [], 10, true);
sd_small = datasim.SystemData(A,B,X_small,Z_small,U_small,data.Phi11,data.Phi12,data.Phi22);
[~, K_small, ~, ~, ~] = sd_small.solve_sdp_on_data();
rho_small = max(abs(eig(A+B*K_small)));
fprintf('  rho_small: %e (元: %e)\n', rho_small, rho_ori);
if ~isempty(history_small) && isfield(history_small, 'rho')
    fprintf('  履歴: rho[0]=%e -> rho[%d]=%e\n', history_small.rho(1), length(history_small.rho)-1, history_small.rho(end));
end

% ============================================
% 結果の比較
% ============================================
fprintf('\n=== 結果の比較 ===\n');
fprintf('元のrho:     %e\n', rho_ori);
fprintf('DIRECT_DGSM_EV:        %e (変化: %+.2e)\n', rho_dgsm, rho_dgsm - rho_ori);
fprintf('IMPLICIT_IDGSM_RHO_LARGE:  %e (変化: %+.2e)\n', rho_large, rho_large - rho_ori);
fprintf('IMPLICIT_IDGSM_RHO_SMALL:  %e (変化: %+.2e)\n', rho_small, rho_small - rho_ori);

% ============================================
% 可視化
% ============================================
fprintf('\n=== 可視化 ===\n');
fprintf('DIRECT_DGSM_EVの結果を表示します...\n');
visualize.plot_data(X, X_dgsm);
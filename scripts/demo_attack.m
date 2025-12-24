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


[X_sdp_adv, Z_sdp_adv, U_sdp_adv] = attack.execute_attack(data, cfg.AttackType.DIRECT_DGSM_EV);

% 差分を計算（グラフ表示用）
dX = X_sdp_adv - X;
dZ = Z_sdp_adv - Z;
dU = U_sdp_adv - U;

% 攻撃後のデータでSystemDataを作成（Phiは元のデータから継承）
sd_sdp_ev = datasim.SystemData(A,B,X_sdp_adv,Z_sdp_adv,U_sdp_adv,data.Phi11,data.Phi12,data.Phi22);
[sol, K, Y, L, diagnostics] = sd_sdp_ev.solve_sdp_on_data();
lambda_sdp_ev = max(abs(eig(A+B*K)));
disp('lambda_sdp_ev');disp(lambda_sdp_ev);

visualize.plot_data(X,X_sdp_adv);
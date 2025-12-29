% U生成→データ生成→SDP→表示の最小デモ
clear; clc; close all;

% 1) 連続の種 or 再現性
% rng(1);

% 2) システム＆重み
[n,m,T] = deal(5,3,cfg.Const.SAMPLE_COUNT);
[A,B] = datasim.make_lti(n,m);

disp('A');disp(A);
disp('B');disp(B);

% 3) 入力とデータ取得
V = make_inputU(m);
[X,Z,U] = datasim.simulate_openloop_stable(A,B,V);


W = Z - A*X - B*U;   % n×T

Phi11 = 1e-1*eye(n);     % 5%マージン
% Phi11 = cfg.Const.SAMPLE_COUNT * eye(n);
Phi12 = zeros(n,T);
Phi22 = -eye(T);


% 4) SDP
data = datasim.SystemData(A,B,X,Z,U,Phi11,Phi12,Phi22);
[sol_ori, K_ori, Y_ori, L_ori, diagnostics_ori] = data.solve_sdp_on_data();
hinf_norm_ori = helper.hinfnorm_AK(A,B,K_ori);
hinf_norm_ori_inv = 1/hinf_norm_ori;


[sol_reg, K_reg, delta_reg, L_reg, diagnostics_reg] = regularization_sdp.solve_sdp(data);
hinf_norm_reg = helper.hinfnorm_AK(A,B,K_reg);
hinf_norm_reg_inv = 1/hinf_norm_reg;

rho_rough = data.rough_quantization_lower_bound();
fprintf('rho_rough: %e\n', rho_rough);
disp("eig(A)");disp(eig(A));

% 5) 表示
delta_ori = sol_ori.delta;
delta_reg = sol_reg.delta;
fprintf('hinf_norm_ori_inv: %e\n', hinf_norm_ori_inv);
fprintf('hinf_norm_reg_inv: %e\n', hinf_norm_reg_inv);
fprintf('delta_ori: %e\n', delta_ori);
fprintf('delta_reg: %e\n', delta_reg);
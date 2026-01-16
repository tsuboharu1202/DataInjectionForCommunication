% U生成→データ生成→SDP→表示の最小デモ
clear; clc; close all;

% 1) システム設定（cfg.Systemから取得）
A = cfg.System.A;
B = cfg.System.B;
[n, m] = cfg.System.getDimensions();
T = cfg.System.getSampleCount();

% disp('A');disp(A);
% disp('B');disp(B);

% 3) 入力とデータ取得
V = make_inputU(m, T);
% [X,Z,U] = datasim.simulate_openloop_stable(A,B,V);
[X,Z] = datasim.simulate_openloop(A,B,V);
U = V;

W = Z - A*X - B*U;   % n×T

Phi11 = 1e-1*eye(n);     % 5%マージン
% Phi11 = T * eye(n);
Phi12 = zeros(n,T);
Phi22 = -eye(T);


% 4) SDP
data = datasim.SystemData(A,B,X,Z,U,Phi11,Phi12,Phi22);
[sol_ori, K_ori, Y_ori, L_ori, diagnostics_ori] = baseline.solve_sdp(data);
hinf_norm_ori = helper.hinfnorm_AK(A,B,K_ori);
hinf_norm_ori_inv = 1/hinf_norm_ori;



% B2 = B + rand(n,m);
rand_mat = rand(n,n);
rand_mat_inv = inv(rand_mat);
B2 = rand_mat*B;
disp("B");disp(B);
disp("B2");disp(B2);
V2 = make_inputU(m, T);
[X2,Z2,U2] = datasim.simulate_openloop_stable(A,B2,V2);
data2 = datasim.SystemData(A,B2,X2,Z2,U2,Phi11,Phi12,Phi22);
[sol_ori2, K_ori2, Y_ori2, L_ori2, diagnostics_ori2] = baseline.solve_sdp(data2);
hinf_norm_ori2 = helper.hinfnorm_AK(A,B2,K_ori2);
hinf_norm_ori2_inv = 1/hinf_norm_ori2;

[sol_reg, K_reg, delta_reg, L_reg, diagnostics_reg] = proposed.solve_sdp(data,1);
hinf_norm_reg = helper.hinfnorm_AK(A,B,K_reg);
hinf_norm_reg_inv = 1/hinf_norm_reg;

rho_rough = data.rough_quantization_lower_bound();
fprintf('rho_rough: %e\n', rho_rough);
disp("eig(A)");disp(eig(A));

disp('Y_ori');disp(value(Y_ori));
disp('Y_ori2');disp(value(Y_ori2));

disp('BK'); disp(B*K_ori);
disp('BK2'); disp(rand_mat_inv*B2*K_ori2);

% 5) 表示
delta_ori = sol_ori.delta;
delta_ori2 = sol_ori2.delta;
delta_reg = sol_reg.delta;
fprintf('hinf_norm_ori_inv1: %e\n', hinf_norm_ori_inv);
fprintf('hinf_norm_ori_inv2: %e\n', hinf_norm_ori2_inv);
fprintf('hinf_norm_reg_inv: %e\n', hinf_norm_reg_inv);
fprintf('delta_ori: %e\n', delta_ori);
fprintf('delta_ori2: %e\n', delta_ori2);
fprintf('delta_reg: %e\n', delta_reg);
disp('K_ori');disp(K_ori);
disp('K_ori2');disp(K_ori2);
disp("K_reg");disp(K_reg);

Gamma_mat = [U;X];
disp('svd(Gamma_mat)');disp(svd(Gamma_mat));
disp('norm(Z,2)');disp(norm(Z,2));
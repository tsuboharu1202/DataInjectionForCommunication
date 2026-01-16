% demo_implicit_regularization.m
% Implicit differentiationによる勾配計算のテスト（正則化版）
%
% 目的：
% 1. deltaの勾配方向に攻撃を加える
% 2. 攻撃前と攻撃後のdeltaを比較する

clear; clc; close all;

fprintf('=== Delta勾配方向への攻撃テスト ===\n\n');

% ============================================
% 1. モデル作成とSDPを解く
% ============================================
fprintf('1. モデル作成とSDPを解く...\n');

% システムサイズ
[n, m, T] = deal(3, 1, cfg.Const.SAMPLE_COUNT);

% システム生成
A = [-0.192, -0.936, -0.814;
    -0.918,  0.729, -0.724;
    -0.412,  0.735, -1.516];
B = [-0.554;
    0.735;
    0.528];

disp('eig A');
disp(eig(A));

% 入力とデータ取得
V = make_inputU(m);
[X, Z, U] = datasim.simulate_openloop_stable(A, B, V);

% Phi設定
Phi11 = 1e-1 * eye(n);
Phi12 = zeros(n, T);
Phi22 = -eye(T);

% 正則化付きSDPを解く
data = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
gamma = 1e3;
[sol, ~, ~, ~, ~] = proposed.solve_sdp(data, gamma);

delta_before = sol.delta;

fprintf('  攻撃前: delta = %.6f\n', delta_before);

% ============================================
% 2. Deltaの勾配を計算
% ============================================
fprintf('\n2. Deltaの勾配を計算...\n');

L_val = sol.L;
Lambda1 = sol.Lambda1;
Lambda2 = sol.Lambda2;
F1 = sol.F1;
Lambda3 = sol.Lambda3;


Gamma = [U; X];
Pi = eye(T) - pinv(Gamma)*Gamma;

dtDelta_dD = implicit_regularization.dtDelta_dD(n,m,T,B,X,Z,U,Pi,gamma,Gamma,L_val,Lambda1,Lambda2,F1,Lambda3);

fprintf('  勾配計算完了\n');

% ============================================
% 3. Deltaの勾配方向に攻撃を加える
% ============================================
fprintf('\n3. Deltaの勾配方向に攻撃を加える...\n');

% 攻撃の強度
attack_strength = 1e-2;

% 攻撃方向：勾配の符号を取る（deltaを小さくする方向 = 負の勾配方向）
attack_direction = -sign(dtDelta_dD);

% ノイズ生成
noise_Z = attack_strength * attack_direction(1:n, :);
noise_X = attack_strength * attack_direction(n+1:2*n, :);
noise_U = attack_strength * attack_direction(2*n+1:end, :);

% 攻撃後のデータ
X_attacked = X + noise_X;
Z_attacked = Z + noise_Z;
U_attacked = U + noise_U;

% ============================================
% 4. 攻撃後のdeltaを計算
% ============================================
fprintf('4. 攻撃後のdeltaを計算...\n');

data_attacked = datasim.SystemData(A, B, X_attacked, Z_attacked, U_attacked, Phi11, Phi12, Phi22);
[sol_attacked, ~, ~, ~, ~] = proposed.solve_sdp(data_attacked, gamma);

delta_after = sol_attacked.delta;
delta_change = delta_after - delta_before;

fprintf('\n=== 結果 ===\n');
fprintf('攻撃前: delta = %.6f\n', delta_before);
fprintf('攻撃後: delta = %.6f\n', delta_after);
fprintf('変化量: Δdelta = %.6e\n', delta_change);
fprintf('相対変化: %.2f%%\n', 100 * delta_change / delta_before);

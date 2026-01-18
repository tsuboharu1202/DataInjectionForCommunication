%% demo_rank_analysis.m
% baselineとproposedの勾配計算でのrank落ち分析
% 一度だけ実行して、どこに自由度があるか確認する

clear; clc;
fprintf('=== Rank落ち分析 ===\n\n');

%% 1. システム生成
% rng(42);  % 再現性のため
n = 3; m = 1; T = 20;

fprintf('システムサイズ: n=%d, m=%d, T=%d\n\n', n, m, T);

% 安定なシステムを生成
A = randn(n)*1.5;
disp('eig(A)');disp(eig(A));
B = randn(n, m);
disp('A');disp(A);
disp('B');disp(B);

% A = cfg.System.A;
% B = cfg.System.B;

% データ生成
x0 = randn(n, 1);
X = zeros(n, T);
Z = zeros(n, T);
U = randn(m, T);
x = x0;
for t = 1:T
    X(:, t) = x;
    x_next = A * x + B * U(:, t);
    Z(:, t) = x_next;
    x = x_next;
end

% Phi行列（IQC用）
Phi11 = 1e-10*eye(n);
Phi12 = zeros(n, T);
Phi22 = -eye(T);

data = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
gamma = 1e3;

%% 2. Proposed手法でのrank分析
fprintf('========================================\n');
fprintf('  PROPOSED手法のrank分析\n');
fprintf('========================================\n');

% SDPを解く
[sol_prop, K_prop, ~, ~, diag_prop] = proposed.solve_sdp(data, gamma);
if diag_prop.problem ~= 0
    error('Proposed SDP failed');
end

% 変数を取得
L = sol_prop.L;
Lambda = sol_prop.Lambda1;
Lambda_P = sol_prop.Lambda2;
Lambda3 = sol_prop.Lambda3;
F1 = sol_prop.F1;

% calc_gradと同じ定義を使う！
Gamma = [U; X];  % [U; X]の順番
Pi = eye(T) - pinv(Gamma)*Gamma;  % 射影行列

% G行列を構築
G1_sol = proposed.gradient.G_grad.G1_merge(n,m,B,T);
G2_sol = proposed.gradient.G_grad.G2_merge(n,m,T,X,Z,U,Pi,gamma,Gamma,L,Lambda,Lambda_P,Lambda3);
G3_sol = proposed.gradient.G_grad.G3_merge(n,m,T,X,Z,U,F1,L,Lambda,B);
G4_sol = proposed.gradient.G_grad.G4_merge(n,m,T);
G5_sol = proposed.gradient.G_grad.G5_merge(n,m,T,X,Lambda_P,L);
G6_sol = proposed.gradient.G_grad.G6_merge(n,m,T);
G7_sol = proposed.gradient.G_grad.G7_merge(n,m,T,X,L);
G8_sol = proposed.gradient.G_grad.G8_merge(n,m,T);

Matrix_H_prop = [
    G1_sol.G1_row_without_Data;
    G2_sol.G2_row_without_Data;
    G3_sol.G3_row_without_Data;
    G4_sol.G4_row_without_Data;
    G5_sol.G5_row_without_Data;
    G6_sol.G6_row_without_Data;
    G7_sol.G7_row_without_Data;
    G8_sol.G8_row_without_Data;
    ];

Matrix_H_prop_full = double(full(Matrix_H_prop));
[n_rows_prop, n_cols_prop] = size(Matrix_H_prop_full);
rank_H_prop = rank(Matrix_H_prop_full);

fprintf('Matrix_H: size=[%d, %d], rank=%d\n', n_rows_prop, n_cols_prop, rank_H_prop);
if rank_H_prop == n_cols_prop
    fprintf('  → 解は一意です (rank = 列数)\n');
else
    fprintf('  → 解は一意ではありません\n');
    fprintf('    自由度: %d (= 列数 %d - rank %d)\n', n_cols_prop - rank_H_prop, n_cols_prop, rank_H_prop);
end

% SVD分析
fprintf('\n--- SVD分析 ---\n');
s_prop = svd(Matrix_H_prop_full);
fprintf('最小特異値: %e\n', min(s_prop));
fprintf('最大特異値: %e\n', max(s_prop));
fprintf('条件数: %e\n', max(s_prop)/min(s_prop));

% 下位特異値
n_show = min(15, length(s_prop));
fprintf('\n下位%d個の特異値:\n', n_show);
for i = 1:n_show
    fprintf('  %2d: %e\n', i, s_prop(end-i+1));
end

% Null space分析
tol_null = 1e-10;
nullvec_prop = null(Matrix_H_prop_full);
fprintf('\nNull space次元: %d (数値的)\n', size(nullvec_prop, 2));

if ~isempty(nullvec_prop)
    % 変数ブロックのサイズ（proposed版）
    n_delta = 1;
    n_L = n * m;
    n_Lambda = (n+1)^2;
    n_Lambda_P = (n+1)^2;
    n_Lambda3 = (n+m)^2;
    
    total_expected = n_delta + n_L + n_Lambda + n_Lambda_P + n_Lambda3;
    fprintf('\n変数の総次元: %d (実際の列数: %d)\n', total_expected, n_cols_prop);
    
    fprintf('\n各変数ブロックでのnull space成分:\n');
    for i = 1:min(3, size(nullvec_prop, 2))
        fprintf('  Null vector %d:\n', i);
        idx = 1;
        fprintf('    d_delta:     norm=%e (idx %d)\n', abs(nullvec_prop(idx, i)), idx);
        idx_L = idx + (1:n_L);
        fprintf('    dL:          norm=%e (idx %d-%d)\n', norm(nullvec_prop(idx_L, i)), idx_L(1), idx_L(end));
        idx_Lambda = idx_L(end) + (1:n_Lambda);
        fprintf('    dLambda:     norm=%e (idx %d-%d)\n', norm(nullvec_prop(idx_Lambda, i)), idx_Lambda(1), idx_Lambda(end));
        idx_Lambda_P = idx_Lambda(end) + (1:n_Lambda_P);
        fprintf('    dLambda_P:   norm=%e (idx %d-%d)\n', norm(nullvec_prop(idx_Lambda_P, i)), idx_Lambda_P(1), idx_Lambda_P(end));
        if idx_Lambda_P(end) < n_cols_prop
            idx_Lambda3 = idx_Lambda_P(end) + (1:n_Lambda3);
            fprintf('    dLambda3:    norm=%e (idx %d-%d)\n', norm(nullvec_prop(idx_Lambda3, i)), idx_Lambda3(1), idx_Lambda3(end));
        end
    end
end

%% 3. Baseline手法でのrank分析
fprintf('\n\n========================================\n');
fprintf('  BASELINE手法のrank分析\n');
fprintf('========================================\n');

% SDPを解く
[sol_base, K_base, Y_base, L_base, diag_base] = baseline.solve_sdp(data);
if diag_base.problem ~= 0
    error('Baseline SDP failed');
end

% 必要な変数を取得
% baseline.solve_sdpと同じGの定義を使う
% G は (3n+m) x (n+T) のサイズ
G = [ eye(n),      Z - B*U;
    zeros(n,n),  -X;
    zeros(n,n),  zeros(n,T);
    zeros(m,n+T) ];  % 最後は m 行！
Phi = [Phi11, Phi12; Phi12', Phi22];


% sol_baseから取得
Y_base = sol_base.Y;
L_base_val = sol_base.L;
alpha = sol_base.alpha;
beta = sol_base.beta;
tDelta = sol_base.delta^2;  % delta = sqrt(tDelta) なので

% F1, F2, F3 は build_lmi_blocks で計算
[F1_base, F2, F3] = baseline.build_lmi_blocks(Y_base, L_base_val, alpha, beta, tDelta, G, Phi, data);

% Dual変数
Lambda1 = sol_base.Lambda1;
Lambda3_base = sol_base.Lambda3;
Lambda_alpha = sol_base.Lambda_alpha;
Lambda_beta = sol_base.Lambda_beta;
Lambda_tDelta = sol_base.Lambda_tDelta;
Lambda_Y = sol_base.Lambda_Y;
Y = Y_base;

% G行列を構築
G1_sol_b = baseline.gradient.G_grad.G1_merge(n,m,B,T);
G2_sol_b = baseline.gradient.G_grad.G2_merge(n,m,B,T);
G3_sol_b = baseline.gradient.G_grad.G3_merge(n,m,T,Lambda1,F2,B,G,Phi);
G4_sol_b = baseline.gradient.G_grad.G4_merge(n,m,T);
G5_sol_b = baseline.gradient.G_grad.G5_merge(n,m,T,B);
G6_sol_b = baseline.gradient.G_grad.G6_merge(n,m,T,Lambda1,alpha,F1_base,F2,B,G,Phi);
G7_sol_b = baseline.gradient.G_grad.G7_merge(n,m,T,Lambda3_base,F3);
G8_sol_b = baseline.gradient.G_grad.G8_merge(n,m,T,alpha,Lambda_alpha);
G9_sol_b = baseline.gradient.G_grad.G9_merge(n,m,T,beta,Lambda_beta);
G10_sol_b = baseline.gradient.G_grad.G10_merge(n,m,T,tDelta,Lambda_tDelta);
G11_sol_b = baseline.gradient.G_grad.G11_merge(n,m,T);
G12_sol_b = baseline.gradient.G_grad.G12_merge(n,m,T);
G13_sol_b = baseline.gradient.G_grad.G13_merge(n,m,T);
G14_sol_b = baseline.gradient.G_grad.G14_merge(n,m,T,Lambda_Y,Y);
G15_sol_b = baseline.gradient.G_grad.G15_merge(n,m,T,Lambda_Y);

Matrix_H_base = [
    G5_sol_b.G5_row_without_Data;
    G1_sol_b.G1_row_without_Data;
    G2_sol_b.G2_row_without_Data;
    G3_sol_b.G3_row_without_Data;
    G4_sol_b.G4_row_without_Data;
    G6_sol_b.G6_row_without_Data;
    G7_sol_b.G7_row_without_Data;
    G8_sol_b.G8_row_without_Data;
    G9_sol_b.G9_row_without_Data;
    G10_sol_b.G10_row_without_Data;
    G11_sol_b.G11_row_without_Data;
    G12_sol_b.G12_row_without_Data;
    G13_sol_b.G13_row_without_Data;
    G14_sol_b.G14_row_without_Data;
    G15_sol_b.G15_row_without_Data
    ];

Matrix_H_base_full = double(full(Matrix_H_base));
[n_rows_base, n_cols_base] = size(Matrix_H_base_full);
rank_H_base = rank(Matrix_H_base_full);

fprintf('Matrix_H: size=[%d, %d], rank=%d\n', n_rows_base, n_cols_base, rank_H_base);
if rank_H_base == n_cols_base
    fprintf('  → 解は一意です (rank = 列数)\n');
else
    fprintf('  → 解は一意ではありません\n');
    fprintf('    自由度: %d (= 列数 %d - rank %d)\n', n_cols_base - rank_H_base, n_cols_base, rank_H_base);
end

% SVD分析
fprintf('\n--- SVD分析 ---\n');
s_base = svd(Matrix_H_base_full);
fprintf('最小特異値: %e\n', min(s_base));
fprintf('最大特異値: %e\n', max(s_base));
fprintf('条件数: %e\n', max(s_base)/min(s_base));

% 下位特異値
n_show = min(15, length(s_base));
fprintf('\n下位%d個の特異値:\n', n_show);
for i = 1:n_show
    fprintf('  %2d: %e\n', i, s_base(end-i+1));
end

% Null space分析
nullvec_base = null(Matrix_H_base_full);
fprintf('\nNull space次元: %d (数値的)\n', size(nullvec_base, 2));

if ~isempty(nullvec_base) && size(nullvec_base, 2) > 0
    % 変数ブロックのサイズ（baseline版）
    n_L_b = n*m;
    n_Y_b = n*n;
    n_alpha_b = 1;
    n_beta_b = 1;
    n_tDelta_b = 1;
    n_Lambda1_b = (3*n+m)*(3*n+m);
    n_Lambda3_b = (n+m)*(n+m);
    n_Lambda_alpha_b = 1;
    n_Lambda_beta_b = 1;
    n_Lambda_tDelta_b = 1;
    n_Lambda_Y_b = n*n;
    
    fprintf('\n各変数ブロックでのnull space成分:\n');
    for i = 1:min(3, size(nullvec_base, 2))
        fprintf('  Null vector %d:\n', i);
        idx = 0;
        idx_L = idx + (1:n_L_b);
        fprintf('    dL:             norm=%e\n', norm(nullvec_base(idx_L, i)));
        idx = idx_L(end);
        idx_Y = idx + (1:n_Y_b);
        fprintf('    dY:             norm=%e\n', norm(nullvec_base(idx_Y, i)));
        idx = idx_Y(end);
        fprintf('    dAlpha:         norm=%e\n', abs(nullvec_base(idx+1, i)));
        idx = idx + 1;
        fprintf('    dBeta:          norm=%e\n', abs(nullvec_base(idx+1, i)));
        idx = idx + 1;
        fprintf('    dtDelta:        norm=%e\n', abs(nullvec_base(idx+1, i)));
        idx = idx + 1;
        idx_Lambda1 = idx + (1:n_Lambda1_b);
        fprintf('    dLambda1:       norm=%e\n', norm(nullvec_base(idx_Lambda1, i)));
        idx = idx_Lambda1(end);
        idx_Lambda3 = idx + (1:n_Lambda3_b);
        fprintf('    dLambda3:       norm=%e\n', norm(nullvec_base(idx_Lambda3, i)));
        idx = idx_Lambda3(end);
        fprintf('    dLambda_alpha:  norm=%e\n', abs(nullvec_base(idx+1, i)));
        idx = idx + 1;
        fprintf('    dLambda_beta:   norm=%e\n', abs(nullvec_base(idx+1, i)));
        idx = idx + 1;
        fprintf('    dLambda_tDelta: norm=%e\n', abs(nullvec_base(idx+1, i)));
        idx = idx + 1;
        idx_Lambda_Y = idx + (1:n_Lambda_Y_b);
        fprintf('    dLambda_Y:      norm=%e\n', norm(nullvec_base(idx_Lambda_Y, i)));
    end
end

%% まとめ
fprintf('\n\n========================================\n');
fprintf('  まとめ\n');
fprintf('========================================\n');
fprintf('Proposed: rank落ち = %d (行: %d, 列: %d, rank: %d)\n', ...
    n_cols_prop - rank_H_prop, n_rows_prop, n_cols_prop, rank_H_prop);
fprintf('Baseline: rank落ち = %d (行: %d, 列: %d, rank: %d)\n', ...
    n_cols_base - rank_H_base, n_rows_base, n_cols_base, rank_H_base);

fprintf('\n注意: rank落ちがあっても、dtDelta部分が一意に決まる場合もある。\n');
fprintf('      null spaceでd_deltaやdtDeltaのノルムが0に近いかを確認すること。\n');


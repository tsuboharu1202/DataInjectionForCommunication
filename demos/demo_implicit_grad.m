% demo_implicit_grad.m
% Implicit differentiationによる勾配計算のテスト
%
% 目的：
% 1. SDPを解く
% 2. その解を用いてimplicit differentiationで勾配を計算
% 3. rankなどの問題を確認

clear; clc; close all;

fprintf('=== Implicit Differentiation 勾配計算テスト ===\n\n');

% ============================================
% 1. モデル作成とSDPを解く（demo_sdp.mを参考）
% ============================================
fprintf('1. モデル作成とSDPを解く...\n');

% システムサイズ（このデモは大きめのシステムでテスト）
n = 4; m = 3;  % 小さめのサイズでテスト
T = 20;  % サンプル数
fprintf('  システムサイズ: n=%d, m=%d, T=%d\n', n, m, T);

% システム生成
[A, B] = datasim.make_lti(n, m);
fprintf('  システム生成完了\n');

% 入力とデータ取得
V = make_inputU(m, T);
[X, Z, U] = datasim.simulate_openloop_stable(A, B, V);
fprintf('  データ生成完了\n');

% Phi設定
Phi11 = 1e-7 * eye(n);
Phi12 = zeros(n, T);
Phi22 = -eye(T);

% SDPを解く
data = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
[sol, K, Y, L, diagnostics] = data.solve_sdp_on_data();

if diagnostics.problem ~= 0
    warning('SDPが正常に解けませんでした。status=%d', diagnostics.problem);
    return;
end

fprintf('  SDP解決完了: rho=%e\n', sol.rho);
fprintf('  alpha=%e, beta=%e, tDelta=%e\n', sol.alpha, sol.beta, sol.delta);

% ============================================
% 2. SDPの解から必要な値を取得
% ============================================
fprintf('\n2. SDPの解から必要な値を取得...\n');

% 解の値を取得
Y_val = sol.Y;
L_val = sol.L;
alpha_val = sol.alpha;
beta_val = sol.beta;
tDelta_val = sol.delta;

% G, Phiを構築（solve_sdpと同じ方法）
G = [eye(n), Z - B*U;
    zeros(n, n), -X;
    zeros(n, n), zeros(n, T);
    zeros(m, n+T)];

Phi = [Phi11, Phi12; Phi12', Phi22];

% F1, F2, F3を構築（解の値で）
[F1, F2, F3] = baseline.build_lmi_blocks(Y_val, L_val, alpha_val, beta_val, tDelta_val, G, Phi, data);

fprintf('  G: %dx%d, Phi: %dx%d\n', size(G,1), size(G,2), size(Phi,1), size(Phi,2));
fprintf('  F1: %dx%d, F2: %dx%d, F3: %dx%d\n', size(F1,1), size(F1,2), size(F2,1), size(F2,2), size(F3,1), size(F3,2));

% ============================================
% 3. ラグランジュ乗数の取得
% ============================================
fprintf('\n3. ラグランジュ乗数の取得...\n');

% solve_sdpからdualを取得
if ~isfield(sol, 'Lambda1') || isempty(sol.Lambda1)
    warning('  solにLambda1が含まれていません。solve_sdpがdualを返すように修正してください。');
    return;
end

Lambda1 = sol.Lambda1;
Lambda3 = sol.Lambda3;
Lambda_alpha = sol.Lambda_alpha;
Lambda_beta = sol.Lambda_beta;
Lambda_tDelta = sol.Lambda_tDelta;
Lambda_Y = sol.Lambda_Y;
fprintf('  Lambda1: %dx%d, Lambda3: %dx%d\n', size(Lambda1,1), size(Lambda1,2), size(Lambda3,1), size(Lambda3,2));
fprintf('  Lambda_alpha=%e, Lambda_beta=%e, Lambda_tDelta=%e\n', Lambda_alpha, Lambda_beta, Lambda_tDelta);

% KKT条件の確認（相補性条件）
% Lambda1 * (F1 - alpha*F2) = 0 (近似的に)
tolerance = 1e-6;
F1_minus_F2 = F1 - alpha_val*F2;
compl_slack_1 = Lambda1 * F1_minus_F2;
compl_slack_3 = Lambda3 * (F3 - tolerance*eye(n+m));
fprintf('  相補性条件チェック:\n');
fprintf('    ||Lambda1*(F1-alpha*F2)|| = %e (0に近いべき)\n', norm(compl_slack_1(:)));
fprintf('    ||Lambda3*(F3-tol*I)|| = %e (0に近いべき)\n', norm(compl_slack_3(:)));

% ============================================
% 4. Implicit differentiationで勾配計算
% ============================================
fprintf('\n4. Implicit differentiationで勾配計算...\n');

try
    % dtDelta_dDを計算
    dtDelta_dD_result = implicit.dtDelta_dD(n, m, T, B, G, Phi, ...
        Lambda1, F1, F2, Lambda3, F3, ...
        alpha_val, Lambda_alpha, beta_val, Lambda_beta, tDelta_val, Lambda_tDelta, Lambda_Y, ...
        Y_val);
    
    fprintf('  勾配計算成功！\n');
    fprintf('  dtDelta_dDのサイズ: %dx%d\n', size(dtDelta_dD_result,1), size(dtDelta_dD_result,2));
    fprintf('  期待されるサイズ: %dx%d (D = [Z'', X'', U'']'')\n', 2*n+m, T);
    
    % サイズチェック
    if isequal(size(dtDelta_dD_result), [2*n+m, T])
        fprintf('  ✓ サイズは正しいです\n');
    else
        warning('  サイズが期待と異なります！');
    end
    
catch ME
    rethrow(ME);  % MATLABの標準エラー表示を使用
end

% ============================================
% 5. Rank確認など（デバッグ用）
% ============================================
fprintf('\n5. 行列のrank確認（デバッグ用）...\n');
fprintf('  [注意] Matrix_Hのrank確認は、dtDelta_dD.m内で直接fprintfしてください。\n');
fprintf('         例: fprintf(''Matrix_H rank: %%d / %%d\\n'', rank(Matrix_H), size(Matrix_H,1));\n');

% 勾配の統計情報（スパース行列の場合はfullに変換）
dtDelta_dD_full = full(dtDelta_dD_result);
fprintf('\n6. 勾配の統計情報:\n');
fprintf('  min: %e, max: %e\n', min(dtDelta_dD_full(:)), max(dtDelta_dD_full(:)));
fprintf('  mean: %e, std: %e\n', mean(dtDelta_dD_full(:)), std(dtDelta_dD_full(:)));
fprintf('  norm: %e\n', norm(dtDelta_dD_full(:)));

% NaN/Infチェック
if any(isnan(dtDelta_dD_full(:))) || any(isinf(dtDelta_dD_full(:)))
    warning('  勾配にNaNまたはInfが含まれています！');
else
    fprintf('  ✓ NaN/Infはありません\n');
end

fprintf('\n=== テスト完了 ===\n');

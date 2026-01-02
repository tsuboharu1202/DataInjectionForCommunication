% demo_implicit_regularization.m
% Implicit differentiationによる勾配計算のテスト（正則化版）
%
% 目的：
% 1. 正則化付きSDPを解く
% 2. その解を用いてimplicit differentiationで勾配を計算
% 3. rankなどの問題を確認

clear; clc; close all;

fprintf('=== Implicit Differentiation 勾配計算テスト（正則化版）===\n\n');

% ============================================
% 1. モデル作成とSDPを解く（demo_sdp.mを参考）
% ============================================
fprintf('1. モデル作成とSDPを解く...\n');

% システムサイズ
[n, m, T] = deal(4, 3, cfg.Const.SAMPLE_COUNT);  % 小さめのサイズでテスト
fprintf('  システムサイズ: n=%d, m=%d, T=%d\n', n, m, T);

% システム生成
[A, B] = datasim.make_lti(n, m);
fprintf('  システム生成完了\n');

% 入力とデータ取得
V = make_inputU(m);
[X, Z, U] = datasim.simulate_openloop_stable(A, B, V);
fprintf('  データ生成完了\n');

% Phi設定
Phi11 = 1e-1 * eye(n);
Phi12 = zeros(n, T);
Phi22 = -eye(T);

% 正則化付きSDPを解く
data = datasim.SystemData(A, B, X, Z, U, Phi11, Phi12, Phi22);
gamma = 1e3;
[sol, ~, ~, ~, ~] = regularization_sdp.solve_sdp(data,gamma);

fprintf('  SDP解決完了: rho=%e\n', sol.rho);
fprintf('  delta=%e\n', sol.delta);

% ============================================
% 2. SDPの解から必要な値を取得
% ============================================
fprintf('\n2. SDPの解から必要な値を取得...\n');

% 解の値を取得
L_val = sol.L;
delta_val = sol.delta;
F1 = sol.F1;

% ラグランジュ乗数の取得
if ~isfield(sol, 'Lambda1') || isempty(sol.Lambda1)
    warning('  solにLambda1が含まれていません。solve_sdpがdualを返すように修正してください。');
    return;
end

Lambda1 = sol.Lambda1;
Lambda2 = sol.Lambda2;
Lambda3 = sol.Lambda3;
fprintf('  Lambda1: %dx%d, Lambda2: %dx%d, Lambda3: %dx%d\n', ...
    size(Lambda1,1), size(Lambda1,2), ...
    size(Lambda2,1), size(Lambda2,2), ...
    size(Lambda3,1), size(Lambda3,2));

% ============================================
% 3. 正則化付きSDPの制約条件の確認
% ============================================
fprintf('\n3. 制約条件の確認...\n');

% const_matを再構築（solve_sdpと同じ方法）
Xm_L_Sym = (X*L_val + (X*L_val)')/2;
const_mat = [Xm_L_Sym, (Z*L_val)', zeros(n,m), (U*L_val)';
    Z*L_val, Xm_L_Sym, delta_val*B, zeros(n,m);
    zeros(m,n), delta_val*B', eye(m), zeros(m,m);
    U*L_val, zeros(m,n), zeros(m,m), eye(m)];

% 相補性条件の確認
tolerance = 1e-8;
compl_slack_1 = Lambda1 * const_mat;
compl_slack_3 = Lambda3 * (Xm_L_Sym - tolerance*eye(n));
fprintf('  相補性条件チェック:\n');
fprintf('    ||Lambda1*const_mat|| = %e (0に近いべき)\n', norm(compl_slack_1(:)));
fprintf('    ||Lambda3*(Xm_L_Sym-tol*I)|| = %e (0に近いべき)\n', norm(compl_slack_3(:)));

% ============================================
% 4. Implicit differentiationで勾配計算
% ============================================
fprintf('\n4. Implicit differentiationで勾配計算...\n');
Gamma = [U; X];
Pi = eye(T) - pinv(Gamma)*Gamma;

try
    % dtDelta_dDを計算
    % 注意: dtDelta_dD.mの実装が不完全な可能性があります
    % 必要な変数（X, Z, U, Pi, gamma, Gamma, L, Lambda, F1）が
    % 関数内で定義されていない場合はエラーが出る可能性があります
    dtDelta_dD_result = implicit_regularization.dtDelta_dD(n,m,T,B,X,Z,U,Pi,gamma,Gamma,L_val,Lambda1,Lambda2,F1);
    
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
    fprintf('  エラーが発生しました: %s\n', ME.message);
    fprintf('  dtDelta_dD.mの実装を確認してください。\n');
    fprintf('  必要な変数（X, Z, U, Pi, gamma, Gamma, L, Lambda, F1）が\n');
    fprintf('  関数内で定義されていない可能性があります。\n');
    rethrow(ME);
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

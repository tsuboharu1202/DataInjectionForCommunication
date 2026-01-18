function test_dF_dL_FD(n, m, T, B, X, Z, U, L_val, delta_val)
%TEST_DF_DL_FD dF_dLの有限差分検証
%   test_dF_dL_FD(n, m, T, B, X, Z, U, L_val, delta_val)
%
%   入力:
%       n, m, T: システムサイズ
%       B: n×m 入力行列
%       X: n×T 状態行列
%       Z: n×T 次状態行列
%       U: m×T 入力行列
%       L_val: T×n 制御パラメータ行列（数値）
%       delta_val: スカラー（量子化パラメータ、数値）
%
%   出力: なし（fprintfで結果を表示）
%
%   作成ファイル一覧:
%       - proposed.gradient.helper.F1_of_L.m
%       - proposed.gradient.helper.test_dF_dL_FD.m (このファイル)
%
%   呼び出し例:
%       % demo_implicit_regularization.m または類似のスクリプトで:
%       % solve_sdp実行後:
%       [sol, ~, ~, ~, ~] = proposed.solve_sdp(data, gamma);
%       L_val = sol.L;
%       delta_val = sol.delta;
%
%       % 一行で呼び出し:
%       proposed.gradient.helper.test_dF_dL_FD(n, m, T, B, X, Z, U, L_val, delta_val);
%
%   または、名前空間を使わずに:
%       test_dF_dL_FD(n, m, T, B, X, Z, U, L_val, delta_val);

fprintf('\n=== dF_dL 有限差分検証 ===\n\n');

% ============================================
% サニティチェック
% ============================================
fprintf('--- サニティチェック ---\n');

% F1_of_Lのサイズ確認
F1_test = proposed.gradient.helper.F1_of_L(n, m, T, B, X, Z, U, L_val, delta_val);
fprintf('F1_of_L サイズ: %d×%d (期待: %d×%d)\n', ...
    size(F1_test,1), size(F1_test,2), 2*n+2*m, 2*n+2*m);
if size(F1_test,1) ~= 2*n+2*m || size(F1_test,2) ~= 2*n+2*m
    error('F1_of_Lのサイズが不正です');
end

% dF_dLのサイズ確認
dF_dL = proposed.gradient.helper.dF_dL(n, m, T, X, Z, U);
fprintf('dF_dL サイズ: %d×%d (期待: %d×%d)\n', ...
    size(dF_dL,1), size(dF_dL,2), (2*n+2*m)^2, T*n);
if size(dF_dL,1) ~= (2*n+2*m)^2 || size(dF_dL,2) ~= T*n
    error('dF_dLのサイズが不正です');
end

fprintf('✓ サイズチェック完了\n\n');

% ============================================
% 有限差分テスト
% ============================================
fprintf('--- 有限差分テスト ---\n');

nTrials = 10;
eps_list = [1e-6, 1e-7, 1e-8];
relerr_all = zeros(nTrials, length(eps_list));

% ベースライン: F1(L_val)
F1_base = proposed.gradient.helper.F1_of_L(n, m, T, B, X, Z, U, L_val, delta_val);
vec_F1_base = F1_base(:);

for trial = 1:nTrials
    % ランダム方向 E (T×n) を生成
    E = randn(T, n);
    vec_E = E(:);
    
    fprintf('Trial %d/%d:\n', trial, nTrials);
    
    for eps_idx = 1:length(eps_list)
        eps = eps_list(eps_idx);
        
        % 有限差分: lhs = (vec(F1(L+epsE)) - vec(F1(L))) / eps
        L_perturbed = L_val + eps * E;
        F1_perturbed = proposed.gradient.helper.F1_of_L(n, m, T, B, X, Z, U, L_perturbed, delta_val);
        vec_F1_perturbed = F1_perturbed(:);
        
        lhs = (vec_F1_perturbed - vec_F1_base) / eps;
        
        % 解析的微分: rhs = dF_dL * vec(E)
        rhs = dF_dL * vec_E;
        
        % 相対誤差
        norm_lhs = norm(lhs);
        relerr = norm(lhs - rhs) / max(1, norm_lhs);
        relerr_all(trial, eps_idx) = relerr;
        
        fprintf('  eps=%.1e: relerr=%.3e\n', eps, relerr);
    end
end

% ============================================
% 結果サマリー
% ============================================
fprintf('\n--- 結果サマリー ---\n');
for eps_idx = 1:length(eps_list)
    eps = eps_list(eps_idx);
    relerr_eps = relerr_all(:, eps_idx);
    fprintf('eps=%.1e: min=%.3e, median=%.3e, max=%.3e\n', ...
        eps, min(relerr_eps), median(relerr_eps), max(relerr_eps));
end

fprintf('\n=== 検証完了 ===\n\n');

end

